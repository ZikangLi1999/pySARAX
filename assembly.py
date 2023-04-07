import copy

from .section import Section


class Assembly:

    def __init__(self, typeName, location) -> None:
        """
        An single assembly

        Input
        -----
        typeName: str, could be ["driver", "control", ...]
        location: str, like "02C08"

        Usage
        -----
        ```python
        pass
        ```
        """
        self.id = -1
        self.type = typeName
        self.location = location
        self.sections = []

        # Auxiliary variables
        self.r = -1
        self.k = -1

    def addSection(self, section, bounds=None) -> None:
        """
        Add a Section into this Assembly

        Input
        -----
        section: class Secion
        bounds: tuple, the upper & lower bound of section

        Caution
        -------
        If the bounds of section are not set, sections must be added strictly by order from bottom to top.
        """
        # Check input
        if type(section) is not Section:
            raise TypeError("Input section is not Section type. section={}".format(section))

        # Modify bounds if given
        if bounds is None:
            pass
        elif type(bounds) is tuple:
            section.bounds = bounds
        else:
            raise TypeError("Input bounds is type {}, which should be tuple".format(type(bounds)))
        
        self.sections.append(section)
    
    def setRefPlane(self, refSec, bias=0.0):
        """
        Set the reference plane of this Assembly

        Input
        -----
        refSec: Section or int, the object or order of the Section locating at reference plane
        bias: float, the distance from the lower bound to the reference plane

        Caution
        -------
        Please call this method after all sections added.
        """
        # Check input
        if type(refSec) is int:
            refSecId = refSec
            refSec = self.sections[refSecId]
        elif type(refSec) is Section:
            # Find the ID of reference Section
            for secId, sec in enumerate(self.sections):
                if sec == refSec:
                    refSecId = secId
        
        # Check heights
        errSecs = []
        for sec in self.sections:
            if sec.height is None:
                errSecs.append(sec.name)
        if len(errSecs) > 0:
            raise RuntimeError("The height of Sections: {} have NOT been set.".format(str(errSecs)))
        
        # Auxiliary function
        sumHeight = lambda secs: sum([sec.height for sec in secs])

        # Determine the coordinates of Section bounds
        for secId, sec in enumerate(self.sections):
            if secId < refSecId:
                lowerBound = - sumHeight(self.sections[secId:refSecId]) - bias
                higherBound = - sumHeight(self.sections[secId+1:refSecId]) - bias
            else:
                lowerBound = sumHeight(self.sections[refSecId:secId]) - bias
                higherBound = sumHeight(self.sections[refSecId:secId+1]) - bias
            sec.bounds = (lowerBound, higherBound)

    def copy(self, typeName, location, deepcopy=False):
        """
        Create a copy of this Assembly

        Input
        -----
        typeName: str, the type name of this copy. Use '@ORIGEN' to represent the original type name
        location: str, the location identifier of this copy
        deepcopy: bool, whether replace copy.copy() with copy.deepcopy()
        """
        # Check input
        if type(typeName) is not str:
            raise TypeError("Input typeName is type {}, which should be str.".format(type(typeName)))
        if type(location) is not str:
            raise TypeError("Input location is type {}, which should be str.".format(type(location)))

        typeName = typeName.replace('@ORIGEN', self.type)

        assemblyCopy = Assembly(typeName=typeName, location=location)
        assemblyCopy.id = -1
        if deepcopy:
            assemblyCopy.sections = copy.deepcopy(self.sections)
        else:
            assemblyCopy.sections = copy.copy(self.sections)

        return assemblyCopy

    def toLAVENDER(self, meshgrid, coolant) -> str:
        """
        Convert the data of Assembly to LAVENDER format (keyword FA_type)

        Input
        -----
        meshgrid: arrayLike, the mesh grid of the Core
        coolant: Section, use to fill the blank under/above the assembly
        """
        slist = ["{:<16}".format("FA_type")]
        slist.append("{:<16d}".format(self.id))

        idlist = []
        # Traverse every mesh grid, and determine which section does it locate in
        for idx in range(len(meshgrid) - 1):
            isUnderAssembly = True
            isAboveAssembly = True
            for sec in self.sections:
                # Determine whether the grid is under/above the current assembly
                if isUnderAssembly and meshgrid[idx] >= sec.bounds[0]:
                    isUnderAssembly = False
                if isAboveAssembly and meshgrid[idx+1] <= sec.bounds[1]:
                    isAboveAssembly = False
                
                # Mesh grid is included in Section
                if sec.bounds[0] <= meshgrid[idx] and sec.bounds[1] >= meshgrid[idx+1]:
                    if len(idlist) == 0:
                        idlist.append("{:d}".format(sec.id))
                        break
                    # The current grid is in the same section as the previous several ones
                    elif int(idlist[-1].split('*')[-1]) == sec.id:
                        if '*' in idlist[-1]:
                            idlist[-1] = "{:d}*{:d}".format(1+int(idlist[-1].split('*')[0]), sec.id)
                            break
                        else:
                            idlist[-1] = "{:d}*{:d}".format(2, sec.id)
                            break
                    else:
                        idlist.append("{:d}".format(sec.id))
                        break
            
            # Mesh grid is under/above the current assembly, fill it with the coolant
            if isUnderAssembly:
                idlist.append("{:d}".format(coolant.id))
            elif isAboveAssembly:
                idlist.append("{:d}".format(coolant.id))
            # If the previous section is also coolant, merge them
            if (isUnderAssembly or isAboveAssembly) and len(idlist) > 1:
                if int(idlist[-2].split('*')[-1]) == coolant.id:
                    if '*' in idlist[-2]:
                        idlist[-2] = "{:d}*{:d}".format(1+int(idlist[-2].split('*')[0]), coolant.id)
                    else:
                        idlist[-2] = "{:d}*{:d}".format(2, coolant.id)
                    idlist.pop()
        
        slist.append(" ".join(idlist))
        return ''.join(slist)

    def check(self) -> bool:
        """
        Check the axial structure of Assembly
        """
        self.sections = sorted(self.sections, key=lambda x: x.bounds[0])
        ifpass = True
        prev = None
        for sec in self.sections:
            if prev is None:
                prev = sec.bounds[1]
                continue
            if prev != sec.bounds[0]:
                # warnings.warn("Section {} in Assembly {} does not match its previous one.".format(sec.name, self.location))
                print("Warning: Section [{}] in Assembly [{}] does not match its previous one.".format(sec.name, self.location))
                ifpass = False
            prev = sec.bounds[1]
        
        return ifpass

    def __str__(self) -> str:
        slist = ["Assembly"]
        slist.append("{:<15}{:d}".format("ID:", self.id))
        slist.append("{:<15}{}".format("type:", self.type))
        slist.append("{:<15}{}".format("location:", self.location))
        slist.append("{:<15}{:d}".format("sections:", len(self.sections)))

        return '\n'.join(slist)


class Assemblies:

    def __init__(self, name) -> None:
        """
        A list of assemblies in the same type

        Input
        -----
        name: str, could be ["driver", "control", "safety", ...] or other identifier
        csvPath: pathLike, The path of "Benchmark Materials CSV Files" directory

        Usage
        -----
        1. List-like: 
        ```python
        for assembly in assemblies:
            assert type(assembly) == Assembly
        ```
        2. Dict-like:
        ```python
        >>> assemblies["02C08"]
        <Assembly object at location "02C08">
        ```
        """
        self.name = name
        self.assemblies = []

        # Auxilary variables for __iter__
        self.locations = None
        self.loc = None
        
    def addAssembly(self, assembly) -> None:
        """
        Add an Assembly into this Assemblies

        Input
        -----
        assembly: class Assembly
        """
        if type(assembly) is not Assembly:
            raise TypeError("Input assembly is not type Assembly. assembly={}".format(assembly))
        
        self.assemblies.append(assembly)

    def __len__(self) -> int:
        return len(self.assemblies)

    def __getitem__(self, key) -> Assembly:
        if type(key) is int:
            return self.assemblies[key]
        if type(key) is str:
            for assembly in self.assemblies:
                if key == assembly.name:
                    return assembly
            # If no assembly's name == key, raise a KeyError
            raise KeyError("No assembly found by key={}".format(key))
        else:
            raise KeyError("Invalid key type {}".format(type(key)))
    
    def __iter__(self):
        self.loc = -1
        # self.locations = tuple(self.assemblies.keys())
        return self

    def __next__(self) -> Assembly:
        self.loc += 1
        # if self.loc >= len(self.locations):
        #     raise StopIteration
        # return self.assemblies[self.locations[self.loc]]
        if self.loc >= len(self.assemblies):
            raise StopIteration
        return self.assemblies[self.loc]
    
    def __str__(self) -> str:
        slist = ["Assemlbies"]
        slist.append("{:<15}{}".format("name:", self.name))
        slist.append("{:<15}{}".format("number:", len(self.assemblies)))

        return '\n'.join(slist)
