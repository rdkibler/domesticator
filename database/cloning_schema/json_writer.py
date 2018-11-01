import json

scheme={
	"parts"	:
		{
		"type 1":
			{
			"left_flank"		:"GCATCGTCTCATCGGTCTCACCCT",
			"left_compatibility"	:"type 8",
			"right_flank"		:"AACGTGAGACCTGAGACGGCAT",
			"right_compatibility"	:"type 2"
			},
		"type 2":
			{
			"left_flank"		:"GCATCGTCTCATCGGTCTCAAACG",
			"left_compatibility"	:"type 1",
			"right_flank"		:"TATGTGAGACCTGAGACGGCAT",
			"right_compatibility"	:"type 3",
			},
		"type 3":
			{
			"left_flank"		:"GCATCGTCTCATCGGTCTCATATG",
			"left_compatibility"	:"type 2",
			"right_flank"		:"GGATCCTGAGACCTGAGACGGCAT",
			"right_compatibility"	:"type 4",
			},
                "type 3a":
                        {
                        "left_flank"            :"GCATCGTCTCATCGGTCTCATATG",
                        "left_compatibility"    :"type 2",
                        "right_flank"	        :"GGTTCTTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 3b",
                        },
                "type 3":
                        {
                        "left_flank"            :"GCATCGTCTCATCGGTCTCATTCT",
                        "left_compatibility"    :"type 3a",
                        "right_flank"	        :"GGATCCTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 4",
                        },
		"type 4":
			{
			"left_flank"	 	:"GCATCGTCTCATCGGTCTCAATCC",
			"left_compatibility"	:"type 3",
			"right_flank"		:"GCTGTGAGACCTGAGACGGCAT",
			"right_compatibility"	:"type 5",
			},
                "type 4a":
                        {
                        "left_flank"            :"GCATCGTCTCATCGGTCTCAATCC",
                        "left_compatibility"    :"type 3",
                        "right_flank"           :"TAACTCGAGTGGCTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 4b",
                        },
                "type 4b":
                        {
                        "left_flank"            :"GCATCGTCTCATCGGTCTCATGGC",
                        "left_compatibility"    :"type 4a",
                        "right_flank"           :"GCTGTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 5",
                        },
                "type 5":
                        {
                        "left_flank"	        :"GCATCGTCTCATCGGTCTCAGCTG",
                        "left_compatibility"    :"type 4",
                        "right_flank"	        :"TACATGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 6",
                        },
                "type 6":
                        {
                        "left_flank"	        :"GCATCGTCTCATCGGTCTCATACA",
                        "left_compatibility"    :"type 5",
                        "right_flank"	        :"GAGTTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 7",
                        },
                "type 7":
                        {
                        "left_flank"	        :"GCATCGTCTCATCGGTCTCAGAGT",
                        "left_compatibility"    :"type 6",
                        "right_flank"	        :"CCGATGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 8",
                        },
                "type 8":
                        {
                        "left_flank"	        :"GCATCGTCTCATCGGTCTCACCGA",
                        "left_compatibility"    :"type 7",
                        "right_flank"	        :"CCCTTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 1",
                        },
                "type 8a":
                        {
                        "left_flank"            :"GCATCGTCTCATCGGTCTCACCGA",
                        "left_compatibility"    :"type 7",
                        "right_flank"           :"CAATTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 8b",
                        },
                "type 8b":
                        {
                        "left_flank"            :"GCATCGTCTCATCGGTCTCACAAT",
                        "left_compatibility"    :"type 8a",
                        "right_flank"           :"CCCTTGAGACCTGAGACGGCAT",
                        "right_compatibility"   :"type 1",
                        }
	},
	"patterns_to_avoid" : [],
	"enzymes_to_avoid" : ["BsaI", "BsmBI", "NotI", "BamHI", "BbsI", "BglII", "EcoRI", "PstI", "SpeI", "XbaI", "XhoI"],
	"type_required"	: [True],
}


with open('MoClo-YTK.json', 'w') as json_file:
	json.dump(scheme, json_file)
