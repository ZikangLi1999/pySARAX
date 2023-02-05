"""
Test model ver1

Set axial parameters by setting bounds directly

Author: LZK
Date: 2023-1-25
"""
import os
import sys
sys.path.append("../../lib")

from pySARAX import Material, Section, Assembly, Core

fuel = Material(name='fuel')
fuel.addElement('U', 1.0, 273.0)

spentFuel = Material(name='spent fuel')
spentFuel.addElement('U', 0.5, 273.0)

extension = Material(name='extension')
extension.addElement('Fe', 2.0, 273.0)

coolant = Material(name='coolant')
coolant.addElement('Na', 1.0, 273.0)


fuelSec = Section(name='fuel', ring=1, pitch=1.0, height=2.0)
fuelSec.appendRod(1.0, fuel)
fuelSec.appendRegion(2.0, coolant)

fuelSecControl = fuelSec.copy(name='fuel in control')

spentFuelSec = Section(name='spent fuel', ring=1, pitch=1.0, height=2.0)
spentFuelSec.appendRod(1.0, spentFuel)
spentFuelSec.appendRegion(2.0, coolant)

extenSec = Section(name='extension', ring=1, pitch=1.0)
extenSec.appendRod(1.2, extension)
extenSec.appendRegion(2.0, coolant)

lowerExtenSecDriver = extenSec.copy(name='lower extension driver')
lowerExtenSecControl = extenSec.copy(name='lower extension control')
upperExtenSecDriver = extenSec.copy(name='upper extension driver')
upperExtenSecControl = extenSec.copy(name='upper extension control')

lowerExtenSecDriver.height = 1.0
lowerExtenSecControl.height = 0.75
upperExtenSecDriver.height = 1.0
upperExtenSecControl.height = 0.75

blanketSec = Section(name='blanket', ring=1, pitch=1.0)
blanketSec.appendRod(1.0, spentFuel)
blanketSec.appendRegion(2.0, coolant)


driver = Assembly(typeName='driver', location='01A01')
driver.addSection(lowerExtenSecDriver)
driver.addSection(fuelSec)
driver.addSection(upperExtenSecDriver)
driver.setRefPlane(fuelSec, bias=1.0)

control = Assembly(typeName='control', location='02A01')
control.addSection(lowerExtenSecControl)
control.addSection(fuelSecControl)
control.addSection(upperExtenSecControl)
control.setRefPlane(1, bias=1.0)

blanket = driver.copy(typeName='blanket', location='02A02')
blanket.sections[1] = spentFuelSec
blanket.sections[1].bounds = (-1.0, 1.0)

core = Core(name='case2', ring=2, pitch=2.0, coolant=extenSec)
core.addRing([driver])
core.addRing([driver, blanket, control, blanket, control, blanket])
# core.complete()
core.specifyID()
core.completed = True
core.meshing(tolerance=0.01)

with open(os.path.join(os.getcwd(), 'TPmate.inp'), 'w', encoding='utf-8') as f:
    tulip = core.toTULIP()
    f.write(tulip)

with open(os.path.join(os.getcwd(), 'lavender.inp'), 'w', encoding='utf-8') as f:
    lavender = core.toLAVENDER()
    f.write(lavender)

core.plotRaial()
core.plotAxial()
