"""
DataBox for LoongSARAX Verification on BER-II Banchmark

Author: LZK
File Log:
2022-12-9  File created
2023-1-23  Rewrite _PERIODIC_TABLE based on periodicTable.PeriodicTable
"""
from .periodicTable import PeriodicTable

"""# Old version: Misused data.
# The Periodic Table using to convert ZAIDS to isotope symbol
periodicTablePath = "C:\SJTUGraduate\Research\Projects\LoongSARAXVerif\code\pySARAX\dataset\各类组件燃料组分\periodicTable.csv"
_PERIODIC_TABLE = pd.read_csv(periodicTablePath)"""

__PROGRAM__ = 'pySARAX'
__AUTHOR__ = "Zikang LI"
__VERSION__ = 1.0

# The Periodic Table using to convert ZAIDS to isotope symbol
_PERIODIC_TABLE = PeriodicTable()

# Group structure of TULIP
_GROUP_STRUCTURE = (1968, 33, 1)

# Temperature of bulk core
_BULK_TEMP = 616.0 # K

# Format data (NOT important)
_INDENT = 4*" "
