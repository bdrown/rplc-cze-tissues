HIT_LISTS=$(wildcard tdreport_dumps/*Hits_List.csv)
HIT_LISTS_PROPS=$(patsubst tdreport_dumps/%Hits_List.csv, tdreport_dumps/%Hits_List_Props.csv, $(HIT_LISTS))

.PHONY : calc_props all

all : calc_props

calc_props : $(HIT_LISTS_PROPS)
	echo "Calculated props"

tdreport_dumps/%Hits_List_Props.csv : tdreport_dumps/%Hits_List.csv calc_protein_props.py
	python calc_protein_props.py $< $@

initialize_python :
	conda activate hubmap_analysis
