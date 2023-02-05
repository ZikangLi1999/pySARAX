"""
Test of the operators *, + of Material

Author: LZK
Date: 2023-1-30
"""
import sys
sys.path.append("C:\\SJTUGraduate\\Research\\Projects\\LoongSARAXVerif\\code\\pySARAX\\lib")

from pySARAX import Material

# Test 1
fuel = Material(name='fuel')
fuel.addNuclide('U235', 1.0, 273.0)

spentFuel = Material(name='spent fuel')
spentFuel.addNuclide('U235', 0.5, 273.0)

result1 = 0.5 * fuel + 1.0 * spentFuel
print(result1)
print(result1.composition)
print(20 * '=')

# Test 2
extension = Material(name='extension')
extension.addElement('Fe', 2.0, 273.0)

coolant = Material(name='coolant')
coolant.addElement('Na', 1.0, 273.0)

result2 = 2.0 * extension + 2.0 * coolant
print(result2)
print(result2.composition)
