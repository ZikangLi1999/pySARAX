import os
import copy
import datetime
import numpy as np

from .assembly import Assembly, Assemblies
from .config import _INDENT, _GROUP_STRUCTURE, __PROGRAM__, __AUTHOR__, __VERSION__


class Core:

    def __init__(self, name, ring, pitch, coolant) -> None:
        """
        A hexagonal core

        Input
        -----
        name: str, the name of this core
        ring: int, the ring number of this core
        pitch: float, the pitch of assembly lattice
        coolant: Section, used to fill the blank at the bottom/top of the Core
        """
        self.name = name
        self.ring = ring
        self.lattice = []
        self.coolant = coolant

        # CONTROL parameters
        self.power = -1.
        self.gammaHeat = False
        self.depletion = False
        self.rodSearch = False
        self.targetKeff = 0.
        self.worth = False
        self.reactivity = False
        self.neutronBalance = False
        self.outputPower = False
        self.outputVTK = True
        self.restart = False
        self.thermal = False
        self.reconstruct = False

        # METHOD parameters
        self.SNOrder = 4
        self.maxInnerIter = 16
        self.maxOuterIter = 500
        self.RMSInner = 5.0E-6
        self.RMSOuter = 1.0E-5
        self.RMSEigen = 1.0E-5
        self.openMPThread = 8
        self.cmacc = True
        self.spaceOrder = 2

        # MATERIAL parameters
        self.lowScatter = False

        # GEOMETRY parametera
        self.meshgrid = np.zeros(shape=0)
        self.boundaryConditions = (0, 0)
        self.pitch = pitch

        # Auxiliary variables
        self.completed = False
        self.assemblies = []
        self.sections = []
        self.materials = []

    def addRing(self, ring) -> None:
        """
        Add a ring of assemblies into this core

        Input
        -----
        ring: list (or Assemblies @TODO), the assemblies in clockwise order, starting at 60 degree (upper-right)
        """
        # Check input
        for assembly in ring:
            assert type(assembly) is Assembly
        
        self.lattice.append(ring)

    def complete(self) -> None:
        """
        Complete Mesh & ID Generation after basic modeling
        """
        if not self.check():
            raise ValueError("Core check not passed. Please check up model again.")
        
        self.specifyID()
        self.completed = True

    def check(self) -> bool:
        """
        Check the structure of this core
        """
        ifpass = True
        # Check whether self.lattice mathes self.ring
        if self.ring != len(self.lattice):
            ifpass = False
            # warnings.warn("Warning: lattice dismatch ring-number of core. Now: {:<6d} Ref: {:<6d}".format(len(self.lattice), self.ring))
            print("Warning: lattice dismatch ring-number of core. Now: {:<6d} Ref: {:<6d}".format(len(self.lattice), self.ring))
        
        # Check the number of assemblies in each ring
        for r, ring in enumerate(self.lattice):
            ref = 6*r if r>0 else 1
            if ref != len(ring):
                ifpass = False
                # warnings.warn("Warning: lattice dismatch hexagonal lattice in ring {:d}. Now {:<6d} Ref: {:<6d}".format(r+1, len(ring), ref))
                print("Warning: lattice dismatch hexagonal lattice in ring [{:d}]. Now {:<6d} Ref: {:<6d}".format(r+1, len(ring), ref))
            for k, assembly in enumerate(ring):
                if not assembly.check():
                    # warnings.warn("Warning: assembly check not passed in ({:d}, {:d})".format(r+1, k+1))
                    print("Warning: assembly check not passed in {}".format((r+1, k+1)))
                    ifpass = False

        return ifpass

    def meshing(self, tolerance=0.1) -> np.ndarray:
        """
        Generate the axial mesh of the whole core

        Input
        -----
        tolerance: float, the maximum height difference to neglect for approximation, ~0.1 cm for SARAX

        Caution
        -------
        The heights of Sections may be modified if there is any Section whose height is close enough to it. 
        The reason is that nodal method calculation may diverge or give wrong result with too-fine mesh grid.
        Therefore, it is suggensted to set the parameter tolerance higher than the minimum distinguishable distance.
        For LoongSARAX, the minimum distinguishable distance is about 0.1 cm, according to Aixin LI in XJTU.
        """
        # Prepare data: collect axial heights of all sections
        heights = set()
        for sec in self.sections:
            heights.add(sec.bounds[0])
            heights.add(sec.bounds[1])
        heights = sorted(list(heights))

        # Approximation: two heights whose distance <= tolerance will be regarded as the same one
        meshgrid = set()
        approxBuffer = set()
        for height in heights:
            # Initialize
            if len(approxBuffer) == 0:
                approxBuffer.add(height)
                continue
            else:
                lower, higher = min(approxBuffer), max(approxBuffer)
                # The current height is close to its adjacent ones, add it into approximation buffer
                if abs(height - lower) <= tolerance and abs(height - higher) <= tolerance:
                    approxBuffer.add(height)
                    continue
                # The current height is far away from its adjacent ones
                # Do approximation on buffer, pop them out to prepare for next approximation
                else:
                    averageHeight = sum(approxBuffer) / len(approxBuffer)
                    meshgrid.add(averageHeight)
                    approxBuffer.clear()
                    approxBuffer.add(height)
        # Do not forget the possible height left in buffer
        if len(approxBuffer) > 0:
            averageHeight = sum(approxBuffer) / len(approxBuffer)
            meshgrid.add(averageHeight)
            approxBuffer.clear()

        # Convert the meshgrid into numpy.ndarray, and sort it
        meshgrid = np.array(list(meshgrid), dtype=np.float)
        meshgrid = np.sort(meshgrid)

        # Approximate the too-close Sections
        for gridHeight in meshgrid:
            for sec in self.sections:
                for i in range(2):
                    if abs(gridHeight - sec.bounds[i]) > 0 and abs(gridHeight - sec.bounds[i]) <= tolerance:
                        bounds = list(sec.bounds)
                        originalBound = bounds[i]
                        bounds[i] = gridHeight
                        sec.bounds = tuple(bounds)
                        print("Warning: The No.{:d} bound {:.4f} of Section ({}) has been approximated to {:.4f}".format(
                            i, originalBound, sec.name, gridHeight))

        self.meshgrid = meshgrid
        self.coolant.bounds = (float(min(self.meshgrid)), float(max(self.meshgrid)))
        return meshgrid

    def specifyID(self) -> None:
        """
        Specify the ID of every Assembly, Section, and Material in this Core
        """
        assemblies = set()
        sections = set()
        materials = set()

        # Traverse the whole Core
        for r, ring in enumerate(self.lattice):
            for k, assembly in enumerate(ring):
                assembly.r = r
                assembly.k = k
                assemblies.add(assembly)

                for section in assembly.sections:
                    if section.eqMethod == 'supercell':
                        sections.add(section.scSection) # Specify ID for supercell section, or core check will NOT pass
                    section.r = r
                    section.k = k
                    sections.add(section)
                    for _, material in section.rod:
                        materials.add(material)
                    for _, material in section.region:
                        materials.add(material)

        sections.add(self.coolant)
        for _, material in self.coolant.region:
            materials.add(material)
        
        # Specify ID
        assemblies = sorted(list(assemblies), key=lambda x: x.location)
        sections = sorted(list(sections), key=lambda x: ' '.join((x.eqMethod, x.name)))
        materials = sorted(list(materials), key=lambda x: x.name)

        for i, assembly in enumerate(assemblies):
            assembly.id = i + 1
        for i, section in enumerate(sections):
            section.id = i + 1
        for i, material in enumerate(materials):
            material.id = i + 1
        
        # Prepare for plot
        self.assemblies = assemblies
        self.sections = sections
        self.materials = materials

    def toTULIP(self) -> str:
        """
        Convert the whole Core to TULIP format
        """
        # Check the integrety of Core
        if not self.completed:
            raise RuntimeError("Core has not completed yet. Please call Core.complete() first.")
        
        # Traverse every sections and materials, and call their methods toTULIP()
        tulipList = []

        # INFO
        tulipList.append("! This file is generated by {}.".format(__PROGRAM__))
        tulipList.append("! CASENAME: {}".format(self.name))
        tulipList.append("! DATE: {}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        tulipList.append("\n")

        # CONTROL
        tulipList.append("! " + 16*"=" + " CONTROL " + 16*"=")
        tulipList.append("CONTROL:")
        tulipList.append("{:<22}{:d}".format("n_mat", len(self.sections)))
        tulipList.append("{:<22}{:d}".format("n_material", len(self.materials)))
        tulipList.append("{:<22}{}".format("leakage_correct", "none"))
        tulipList.append("{:<22}{:<8d}{:<8d}{:<8d}".format("group_info", *_GROUP_STRUCTURE))

        geom_kind = []
        for sec in self.sections:
            if sec.eqMethod == "homo":
                gk = "1"
            elif sec.eqMethod == "1-D":
                gk = "2"
            elif sec.eqMethod == "supercell":
                gk = "3"
            else:
                raise ValueError("Incorrect equivalence method {} in section {}".format(sec.eqMethod, sec.name))
            # Use "xx*yy" form if the last one is same as the current one
            if len(geom_kind) > 0 and geom_kind[-1][-1] == gk:
                if '*' not in geom_kind[-1]:
                    geom_kind[-1] = "2*{}".format(gk)
                else:
                    geom_kind[-1] = str(1 + int(geom_kind[-1][:-2])) + '*' + gk
            else:
                geom_kind.append(gk)
        geom_kind = ' '.join(geom_kind)
        tulipList.append("{:<22}{}".format("geom_kind", geom_kind))
        tulipList.append("\n")
        
        # GEOMETRY
        tulipList.append("! " + 16*"=" + " GEOMETRY " + 16*"=")
        tulipList.append("GEOMETRY:")
        for section in self.sections:
            tulipList.append(section.toTULIP())
            # Supercell Equivalence Method
            if section.eqMethod == "supercell":
                # Locate the adjacent assembly and append it as "sc_rod" and "sc_region"
                tulipList.append(section.scSection.toTULIP(sc=True))
        
        # MATERIAL
        tulipList.append("! " + 16*"=" + " MATERIAL " + 16*"=")
        tulipList.append("MATERIAL:")
        for material in self.materials:
            tulipList.append(material.toTULIP())
        
        return "\n".join(tulipList)

    def toLAVENDER(self) -> str:
        """
        Convert core data to LAVENDER format
        """
        # Check the integrety of Core
        if not self.completed:
            raise RuntimeError("Core has not completed yet. Please call Core.complete() first.")

        slist = []

        # INFO
        slist.append("! This file is generated by {}.".format(__PROGRAM__))
        slist.append("! CASENAME: {}".format(self.name))
        slist.append("! DATE: {}".format(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
        slist.append("\n")

        # CASENAME
        slist.append("! " + 36*"=" + " CASENAME " + 36*"=")
        slist.append("CASENAME: {}".format(self.name))
        slist.append("\n")

        # CONTROL
        bool2tf = lambda x: "T" if x else "F"
        slist.append("! " + 36*"=" + " CONTROL " + 36*"=")
        slist.append("CONTROL:")
        slist.append("{:<16}{:<16}{:<16}".format("!", "power(Wth)", "(n,gamma)heating"))
        slist.append("{:<16}{:<16.4E}{:<16}".format("steady", self.power, bool2tf(self.gammaHeat)))
        slist.append("{:<16}{:<16}".format("!", "on-off"))
        slist.append("{:<16}{:<16}".format("depletion", bool2tf(self.depletion)))
        slist.append("{:<16}{:<16}{:<16}".format("!", "on-off", "target keff"))
        slist.append("{:<16}{:<16}{:<16.5f}".format("rod_search", bool2tf(self.rodSearch), self.targetKeff))
        slist.append("{:<16}{:<16}".format("!", "on-off"))
        slist.append("{:<16}{:<16}".format("worth", bool2tf(self.worth)))
        slist.append("{:<16}{:<16}{:<16}".format("!", "on-off", "neutron balance"))
        slist.append("{:<16}{:<16}{:<16}".format("reactivity", bool2tf(self.reactivity), bool2tf(self.neutronBalance)))
        slist.append("{:<16}{:<16}{:<16}{:<16}".format("!", "power", "VTK", "restart"))
        slist.append("{:<16}{:<16}{:<16}{:<16}".format("output", bool2tf(self.outputPower), bool2tf(self.outputVTK), bool2tf(self.restart)))
        slist.append("{:<16}{:<16}".format("!", "on-off"))
        slist.append("{:<16}{:<16}".format("th", bool2tf(self.thermal)))
        slist.append("{:<16}{:<16}".format("!", "on-off"))
        slist.append("{:<16}{:<16}".format("reconstruct", bool2tf(self.reconstruct)))
        slist.append("\n")

        # METHOD
        slist.append("! " + 36*"=" + " METHOD " + 36*"=")
        slist.append("METHOD:")
        slist.append("{:<16}{:<16}".format("!", "even integer"))
        slist.append("{:<16}{:<16d}".format("sn_order", self.SNOrder))
        slist.append("{:<16}{:<16}{:<16}{:<16}{:<16}{:<16}".format("!", "max_inner", "max_outer", "rms_inner", "rms_outer", "rms_eigen"))
        slist.append("{:<16}{:<16d}{:<16d}{:<16.1E}{:<16.1E}{:<16.1E}".format("error_eigen", self.maxInnerIter, self.maxOuterIter, self.RMSInner, self.RMSOuter, self.RMSEigen))
        slist.append("{:<16}{:<16}".format("!", "nThread"))
        slist.append("{:<16}{:<16}".format("openmp", self.openMPThread))
        slist.append("{:<16}{:<16}".format("!", "on-off"))
        slist.append("{:<16}{:<16}".format("cmacc", bool2tf(self.cmacc)))
        slist.append("{:<16}{:<16}".format("!", "spcae expansion order"))
        slist.append("{:<16}{:<16d}".format("space_order", self.spaceOrder))
        slist.append("\n")

        # MATERIAL (Actually SECTION)
        slist.append("! " + 36*"=" + " MATERIAL " + 36*"=")
        slist.append("MATERIAL:")
        slist.append("{:<16}{:<16}".format("!", "on-off"))
        slist.append("{:<16}{:<16}".format("low_scat_order", bool2tf(self.lowScatter)))
        slist.append("{:<16}{:<16}{:<16}{:<16}{:<16}".format("!", "mat id", "xs file", "temp id", "mat type"))
        for mat in sorted(self.sections, key=lambda x: x.id):
            slist.append("{:<16}{:<16d}{:<16}{:<16d}".format("mat_file", mat.id, "MAT{:d}".format(mat.id), 1))
        slist.append("\n")

        # GEOMETRY
        slist.append("! " + 36*"=" + " GEOMETRY " + 36*"=")
        slist.append("GEOMETRY:")

        # GEOMETRY: LAYER
        sLayer = []
        for grid in np.diff(self.meshgrid):
            sLayer.append("{:.4f}".format(grid))
        slist.append("{:<16}{:<16}{:<16}".format("!", "num_layer", "axial height"))
        slist.append("{:<16}{:<16d}{:<16}".format("layer", len(sLayer), " ".join(sLayer)))

        # GEOMETRY: BC_AXIAL
        slist.append("{:<16}{:<16}{:<16}".format("!", "bc bottom", "bc top"))
        slist.append("{:<16}{:<16d}{:<16d}".format("bc_axial", *self.boundaryConditions))

        # GEOMETRY: FA_TYPE
        slist.append("\n")
        slist.append("{:<16}{:<16}{:<32}{:<32}".format("!", "FA id", "assembly type", "assembly location"))
        for assembly in sorted(self.assemblies, key=lambda x: x.id):
            slist.append("{:<16}{:<16d}{:<32}{:<32}".format("!", assembly.id, assembly.type, assembly.location))
        slist.append("!" + 100*".")
        slist.append("{:<16}{:<16}{:<16}".format("!", "FA id", "layer mat id"))
        for assembly in sorted(self.assemblies, key=lambda x: x.id):
            slist.append(assembly.toLAVENDER(self.meshgrid, self.coolant))
        slist.append("\n")
        
        # GEOMETRY: HEX_DIM
        slist.append("{:<16}{:<16}{:<16}{:<16}".format("!", "degree", "ring", "pitch [cm]"))
        slist.append("{:<16}{:<16d}{:<16d}{:<16.4f}".format("hex_dim", 360, self.ring, self.pitch))

        # GEOMETRY: HEX_CONF
        slist.append("{:<16}{:<16}".format("!", "assembly lattice"))
        slist.append("hex_conf")

        # Traverse the core row-by-row
        #     l: row    r: ring    k: clockwise order
        for l in range(2*self.ring - 1):
            dis2centerLine = abs(self.ring - 1 - l)
            rowlist = [_INDENT, "{:>2d}".format(l+1), _INDENT, (dis2centerLine)*' ']

            # Go from outer ring to inner ring (120 ~ 240 degree)
            for j in range(self.ring - 1 - dis2centerLine):
                r = (self.ring - 1) - j
                k = 5 * r + j - l
                rowlist.append("{:>4}".format(self.lattice[r][k].id))
            
            # Go through the the inneset ring in current row (60 ~ 120 & 240 ~ 300 degree)
            r = dis2centerLine
            kCenterLine = 4 * r if r > 0 else 1
            for j in range(dis2centerLine):
                k = kCenterLine + (r + j) * int(np.sign(self.ring - 1 - l))
                rowlist.append("{:>4}".format(self.lattice[r][k].id))
            
            # The assembly in line of 60 & 300 degree
            k = 2 * dis2centerLine if self.ring - 1 - l < 0 else 0
            rowlist.append("{:>4}".format(self.lattice[r][k].id))

            # Go from the inner ring to outer ring (300 ~ 60 degree)
            for j in range(self.ring - 1 - dis2centerLine):
                r = 1 + j + dis2centerLine
                k = l + r + 1 - self.ring
                rowlist.append("{:>4}".format(self.lattice[r][k].id))
            
            slist.append(''.join(rowlist))
        slist.append("\n")

        # END
        slist.append("! " + 36*"=" + " END " + 36*"=")
        slist.append("END:")

        return '\n'.join(slist)

    def plotRaial(self, pitch=1.0, savePath=None):
        """
        Plot the radial structure of Core

        Input
        -----
        pitch: float, the pitch of core lattice
        savePath: pathLike, the path to save the image, ending with ".png", ".jpg", ".svg" or other supported format
        """
        # Check the completeness of Core
        if not self.completed:
            raise RuntimeError("Core has not completed yet. Please call Core.complete() first.")
        
        # Import matplotlib
        try:
            import matplotlib.pyplot as plt
            from matplotlib import patches
        except ImportError:
            # warnings.warn("Module matplotlib.pyplot could not be imported.")
            print("Warning: Module matplotlib.pyplot could not be imported.")
            return

        # Prepare data for plot
        hBound = pitch * (self.ring+1)
        vBound = pitch * (self.ring+1) * np.cos(np.pi/6)

        # Traverse every assembly and give color by their types
        xcoordinates = []
        ycoordinates = []
        ccoordinates = []
        assemblyTypes = set()
        for r in range(self.ring):
            numInRing = 6 * r if r > 0 else 1
            for k in range(numInRing):
                assembly = self.lattice[r][k]
                assemblyTypes.add(assembly.type)
                x, y = self.rk2xy(r, k, pitch=pitch)
                xcoordinates.append(x)
                ycoordinates.append(y)
                ccoordinates.append(assembly.type)
        
        # Convert type name to integer number
        assemblyTypes = list(assemblyTypes)
        for i, typeName in enumerate(assemblyTypes):
            for j in range(len(ccoordinates)):
                if ccoordinates[j] == typeName:
                    ccoordinates[j] = i + 1 # Color is white if C==0

        # The coordinates of the coolant outside of assemblies lattice (not exact)
        radius = pitch * self.ring
        coolantCirc = patches.Circle(xy=(0., 0.), radius=radius, fill=False)

        # Plot
        fig, axes = plt.subplots()
        axes.add_artist(coolantCirc)
        plt.hexbin(xcoordinates, ycoordinates, ccoordinates, gridsize=(2*self.ring-1-self.ring%2, self.ring-1))
        # plt.scatter(xcoordinates, ycoordinates, alpha=0.3)
        plt.xlim((-hBound, hBound))
        plt.ylim((-vBound, vBound))

        # Save
        if savePath:
            diretory, file = os.path.split(savePath)
            if not os.path.exists(diretory):
                # warnings.warn("The directory [{}] does not exist.".format(diretory))
                print("Warning: The directory [{}] does not exist.".format(diretory))
            else:
                plt.savefig(savePath)
        else:
            plt.show()

    def plotAxial(self, savePath=None):
        """
        Plot the axial structure of all types of Assembly (NOT Core)

        Input
        -----
        savePath: pathLike, the path to save the image, ending with ".png", ".jpg", ".svg" or other supported format
        """
        # Check the completeness of Core
        if not self.completed:
            raise RuntimeError("Core has not completed yet. Please call Core.complete() first.")
        
        # Import Matplotlib
        try:
            import matplotlib.pyplot as plt
            import matplotlib.colors as mcolors
        except ImportError:
            # warnings.warn("Module matplotlib.pyplot could not be imported.")
            print("Warning: Module matplotlib.pyplot could not be imported.")
            return
        
        # Create a new canvas, to prevent paint in the radial one
        plt.figure()
        plt.xticks(rotation=90)
        
        # Prepare data: one assembly represents one type
        assemblies = {}
        for assembly in self.assemblies:
            if assembly.type not in assemblies:
                assemblies[assembly.type] = assembly
        
        # Prepare data: count the total number of sections in representative assemblies
        sections = set()
        colors = list(mcolors.CSS4_COLORS.keys()) # 148 colors
        cmap = {}
        for assembly in assemblies.values():
            for sec in assembly.sections:
                sections.add(sec)
        for idx, sec in enumerate(sections):
            idx = idx % 148 # Prevent the idx >= 148, although it has a low probabliity
            cmap[sec] = colors[idx]

        # Plot: x-axis is assembly type, y-axis is the height of axial structure
        for typeName, assembly in assemblies.items():
            for sec in assembly.sections:
                lower, higher = sec.bounds
                if lower * higher >= 0:
                    if lower >= 0 and higher >= 0:
                        bottom, top = lower, higher
                    elif lower < 0 and higher <= 0:
                        bottom, top = higher, lower
                    plt.bar(x=[typeName], height=[top - bottom], bottom=[bottom], color=[cmap[sec]])
                else:
                    # The signs of two bounds are different, split this section into two
                    bottom, top = 0., higher
                    plt.bar(x=[typeName], height=[top], bottom=[bottom], color=[cmap[sec]])
                    bottom, top = lower, 0.
                    plt.bar(x=[typeName], height=[bottom], bottom=[top], color=[cmap[sec]])
        
        # Plot the reference plain
        plt.axhline(y=0.)
        # Plot the mesh grid
        for grid in self.meshgrid:
            plt.axhline(y=grid, alpha=0.3)

        # Save
        if savePath:
            diretory, file = os.path.split(savePath)
            if not os.path.exists(diretory):
                # warnings.warn("The directory [{}] does not exist.".format(diretory))
                print("Warning: The directory [{}] does not exist.".format(diretory))
            else:
                plt.savefig(savePath)
        else:
            plt.show()

    def rk2xy(self, r, k, pitch=1.0) -> tuple:
        """
        Convert the (ring,clock) coordinate to (x,y) coordinate

        Input
        -----
        r: int, the ring coordinate
        k: int, the clock coordinate
        pitch: float, the pitch of core lattice
        """
        # Check input
        if type(r) is not int or type(k) is not int:
            raise TypeError("Input coordinate is not type int.")
        
        # Handle the innest ring (single assembly)
        if r == 0:
            assert k == 0
            return (0., 0.)

        # Determine which sector is (r,k) in. 0: A, 1: B, ...
        sector = k // r
        # Determine which seat is (r,k) in this sector.
        seat = k % r

        # This first seat in sector (clockwise)
        degree = (1 - sector) * np.pi / 3
        radius = r * pitch
        x = radius * np.cos(degree)
        y = radius * np.sin(degree)

        # Move to the target seat
        degree -= 2 * np.pi / 3
        distance = seat * pitch
        x += distance * np.cos(degree)
        y += distance * np.sin(degree)

        return (x, y)

    def copy(self, name, deepcopy=False):
        """
        Create a copy of this Core

        Input
        -----
        deepcopy: bool, whether replace copy.copy() with copy.deepcopy()
        """
        coreCopy = Core(name=name, ring=self.ring)
        if deepcopy:
            coreCopy.lattice = copy.deepcopy(self.lattice)
        else:
            coreCopy.lattice = copy.copy(self.lattice)
        coreCopy.complete()
        return coreCopy
