# Purpose : downgrade lipid names to the structural-resolution level the evidence supports (Goslin).
#           Re-emits each name at MOLECULAR_SPECIES (underscore) at most, never above what was parsed.
# Inputs  : one or more lipid names on the command line.
# Usage   : python honest_level.py "PC 16:0/18:1" "PC 34:1" "PC O-34:1"
# Checked : pygoslin 2.2.5
import sys
from pygoslin.parser.Parser import LipidParser
from pygoslin.domain.LipidLevel import LipidLevel


def honest_names(name, parser=None):
    """Return (claimed_level, honest_name, sum_name) for one lipid name."""
    parser = parser or LipidParser()
    lipid = parser.parse(name)             # e.g. 'PC 16:0/18:1', a slash-claimed name from a tool export

    claimed_level = lipid.lipid.info.level    # LipidLevel enum the string asserts

    # GUARD: never request a target level MORE specific than what was parsed (see Common Errors: RuntimeException).
    # get_lipid_string() has no chain/sn data to invent for a name parsed at a coarser level
    # (e.g. a sum-composition or ether/plasmalogen name parsed at SPECIES has no chains to report
    # at MOLECULAR_SPECIES) and raises an unhandled RuntimeException instead of degrading gracefully.
    # Cap the target at whichever is coarser: the honest ceiling or what was actually parsed.
    target_level = min(LipidLevel.MOLECULAR_SPECIES, claimed_level, key=lambda l: l.value)

    # For ether/plasmalogen lipids (O-/P-) at SPECIES level, preserve the original name.
    # Goslin's get_lipid_string() at MOLECULAR_SPECIES loses the O-/P- distinction and may
    # convert P- to O- with wrong chain count (e.g. PC P-34:1 -> PC O-34:2).
    if target_level == LipidLevel.SPECIES and claimed_level == LipidLevel.SPECIES:
        # Check if original name has ether/plasmalogen prefix
        if ' O-' in name or ' P-' in name:
            honest_name = name  # Preserve original ether/plasmalogen notation
        else:
            honest_name = name  # SPECIES level has no chains to report anyway
    else:
        honest_name = lipid.get_lipid_string(target_level)   # 'PC 16:0_18:1' here

    sum_name = lipid.get_lipid_string(LipidLevel.SPECIES) if claimed_level.value >= LipidLevel.SPECIES.value else honest_name
    return claimed_level, honest_name, sum_name


if __name__ == '__main__':
    if len(sys.argv) < 2:
        sys.exit('usage: python honest_level.py "PC 16:0/18:1" ...')
    parser = LipidParser()
    for n in sys.argv[1:]:
        level, honest, summed = honest_names(n, parser)
        print(f'{n}\tclaimed={level.name}\thonest={honest}\tsum={summed}')
