import RNA
from inspect import signature

seq = "uaauacgacucacuauaggggaauugugagcggauaacaauuccccucuagaaauaauuuuguuuaacuuuaagaaggagauauacauaugggguccucucaucaucaccaccaucauucgagcggcgaaaaucuguacuuucagggcucaccggaauucuuggagaaccugcgccgccaucucgaccgucugcgcaaacacauccguaaauuagaagaaauccuggaggaaauugacgaucagauggugaaagaagcaaucgaucuuucgcggcgcagcgucgagauuguugaggaagugauugagacguuugaacgugcuggugacgcgaguccgaccaaacuggaugaauuagagaagauucuccaacgugcggaggaaaucauugaucgcgccgauaaacuguuggaauaucgucgccgguaacucgagcaccaccaccaccaccacugagauccggcugcuaacaaagcccgaaaggaagcugaguuggcugcugccaccgcugacaauaacuagcauaaccccuuggggccucuaaacgggucuugagggguuuuuug"

def mfe_window_callback(start,end,structure,energy,data=None):
	data.append({'structure':structure,'start':start,'end':end,'energy':energy})


fc = RNA.fold_compound(seq, None, RNA.OPTION_MFE | RNA.OPTION_WINDOW)
data = []
mfe = fc.mfe_window_cb(mfe_window_callback, data)

max_energy = -5.0

for datum in data:
	if datum['energy'] < max_energy:
		print(datum['start'],datum['end'], datum['end']-datum['start'], datum['energy'])
print(mfe)