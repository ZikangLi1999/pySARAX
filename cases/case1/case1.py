"""
Test case for pySARAX

Author: LZK

Target:
- Test the function of generating various assemblies by ONE template assembly
- Test the function of supercell

File log:
2022-12-15  File created
"""
import sys
sys.path.append("C:\\SJTUGraduate\\Research\\Projects\\LoongSARAXVerif\\code\\pySARAX\\lib")

import os
from pySARAX import Material, Section, Assembly, Assemblies, Core

cwd = "C:\\SJTUGraduate\\Research\\Projects\\LoongSARAXVerif\\code\\pySARAX\\cases\\case1"

# Material
materials = []

uo2 = Material(name="UO2")
uo2.addNuclide(nuclide="U238", density=0.25, temperature=616.0)
uo2.addNuclide(nuclide="O16", density=0.5, temperature=616.0)
materials.append(uo2)

puo2 = Material(name="PuO2")
puo2.addNuclide(nuclide="Pu239", density=0.25, temperature=616.0)
puo2.addNuclide(nuclide="O16", density=0.5, temperature=616.0)
materials.append(puo2)

water = Material(name="water")
water.addNuclide(nuclide="H1", density=0.2, temperature=400.0)
water.addNuclide(nuclide="O16", density=0.1, temperature=400.0)
materials.append(water)

b4c = Material(name="B4C")
b4c.addNuclide(nuclide="B10", density=0.4, temperature=616.0)
b4c.addNuclide(nuclide="C12", density=0.1, temperature=616.0)
materials.append(b4c)

zircaloy = Material(name="zircaloy")
zircaloy.addNuclide(nuclide="Zr91", density=0.9, temperature=550.0)
materials.append(zircaloy)

ss304 = Material(name="SS304")
ss304.addNuclide(nuclide="Fe56", density=0.56, temperature=515.0)
ss304.addNuclide(nuclide="C12", density=0.04, temperature=515.0)


# Section
uoxSlug = Section(name="UOX slug", ring=6, pitch=0.5665, eqMethod="1-D")
uoxSlug.appendRod(size=0.3302, material=uo2)
uoxSlug.appendRod(size=0.3810, material=zircaloy)
uoxSlug.appendRegion(size=0.4591, material=water)
uoxSlug.appendRegion(size=5.8166, material=ss304)
uoxSlug.appendRegion(size=5.8929, material=water)
uoxSlug.bounds = (-1.0, 1.0)

puoxSlug = Section(name="PuOX slug", ring=6, pitch=0.5665, eqMethod="1-D")
puoxSlug.appendRod(size=0.3302, material=puo2)
puoxSlug.appendRod(size=0.3810, material=zircaloy)
puoxSlug.appendRegion(size=0.4591, material=water)
puoxSlug.appendRegion(size=5.8166, material=ss304)
puoxSlug.appendRegion(size=5.8929, material=water)
puoxSlug.bounds = (-1.0, 1.0)

poisonSlug = Section(name="poison slug", ring=3, pitch=0.5665, eqMethod="supercell")
poisonSlug.appendRod(size=0.3302, material=b4c)
poisonSlug.appendRod(size=0.3810, material=zircaloy)
poisonSlug.appendRegion(size=0.4591, material=water)
poisonSlug.appendRegion(size=5.8166, material=ss304)
poisonSlug.appendRegion(size=5.8929, material=water)
poisonSlug.bounds = (-1.0, 1.0)
poisonSlug.scSection = uoxSlug

shieldBlockAbove = Section(name="shield block", ring=0, pitch=0.)
shieldBlockAbove.appendRegion(size=5.8166, material=ss304)
shieldBlockAbove.appendRegion(size=5.8929, material=water)
shieldBlockAbove.bounds = (1.0, 2.0)

shieldBlockBelow = Section(name="shield block", ring=0, pitch=0.)
shieldBlockBelow.appendRegion(size=5.8166, material=ss304)
shieldBlockBelow.appendRegion(size=5.8929, material=water)
shieldBlockBelow.bounds = (-2.0, -1.0)

dummySlug = Section(name="dummy slug", ring=0, pitch=0.)
dummySlug.appendRegion(size=5.8166, material=ss304)
dummySlug.appendRegion(size=5.8929, material=water)
dummySlug.bounds = (-1.0, 1.0)

coolantSection = Section(name="coolant", ring=0, pitch=0.)
coolantSection.appendRegion(size=5.8929, material=water)

# Assembly
driverTemplate = Assembly(typeName="driver", location="template")
driverTemplate.addSection(section=shieldBlockBelow)
driverTemplate.addSection(section=shieldBlockAbove)

control = Assembly(typeName="control", location="01A01")
control.addSection(shieldBlockBelow)
control.addSection(poisonSlug)
control.addSection(shieldBlockAbove)

dummy = Assembly(typeName="dummy", location="03A00")
dummy.addSection(shieldBlockBelow)
dummy.addSection(dummySlug)
dummy.addSection(shieldBlockAbove)

# Core
core = Core(name="case1", ring=3, pitch=5.8929, coolant=coolantSection)

ring1 = Assemblies(name="ring1")
ring1.addAssembly(control)

ring2 = Assemblies(name="ring2")
for k in range(6):
    slug = uoxSlug if k % 2 == 0 else puoxSlug
    driverK = driverTemplate.copy()
    driverK.name = "driver{:d}".format(k)
    driverK.location = "02A{:0>2d}".format(k)
    driverK.addSection(slug)
    ring2.addAssembly(driverK)

ring3 = Assemblies(name="ring3")
for k in range(12):
    dummyK = dummy.copy()
    dummyK.name = "dummy{:d}".format(k)
    dummyK.location = "03A{:0>2d}".format(k)
    ring3.addAssembly(dummyK)

core.lattice = [ring1, ring2, ring3]
core.complete()
core.meshing(tolerance=0.1000)

core.plotRaial(savePath=os.path.join(cwd, "radial.png"))
core.plotAxial(savePath=os.path.join(cwd, "axial.png"))

tulipPath = os.path.join(cwd, "TPmate-case1.inp")
lavenderPath = os.path.join(cwd, "lavender-case1.inp")

with open(tulipPath, 'w', encoding='utf-8') as tulip:
    tulip.write(core.toTULIP())

with open(lavenderPath, 'w', encoding='utf-8') as lavender:
    lavender.write(core.toLAVENDER())
