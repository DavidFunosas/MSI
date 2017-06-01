#!/usr/bin/python

lipid_parent = {
 	'Cholesterol':'Sterols',
	'ERG':'Sterols',
	'DLPA':'PA',
	'DMPA':'PA',
	'DPPA':'PA',
	'DSPA':'PA',
	'POPA':'PA',
	'PLPA':'PA',
	'SOPA':'PA',
	'SLPA':'PA',
	'DYPA':'PA',
	'YOPA':'PA',
	'DOPA':'PA',
	'DGPA':'PA',
	'DEPA':'PA',
	'DNPA':'PA',
	'DDPC':'PC',
	'DCPC':'PC',
	'DLPC':'PC',
	'DMPC':'PC',
	'DPPC':'PC',
	'DSPC':'PC',
	'POPC':'PC',
	'PLPC':'PC',
	'SOPC':'PC',
	'SLPC':'PC',
	'DYPC':'PC',
	'YOPC':'PC',
	'DOPC':'PC',
	'DUPC':'PC',
	'DGPC':'PC',
	'DEPC':'PC',
	'DNPC':'PC',
	'DLPE':'PE',
	'DMPE':'PE',
	'DPPE':'PE',
	'DSPE':'PE',
	'PYPE':'PE',
	'POPE':'PE',
	'PLPE':'PE',
	'SOPE':'PE',
	'SLPE':'PE',
	'DYPE':'PE',
	'YOPE':'PE',
	'OYPE':'PE',
	'DOPE':'PE',
	'DGPE':'PE',
	'DEPE':'PE',
	'DNPE':'PE',
	'DLPG':'PG',
	'DMPG':'PG',
	'DPPG':'PG',
	'DSPG':'PG',
	'PYPG':'PG',
	'POPG':'PG',
	'PLPG':'PG',
	'SOPG':'PG',
	'SLPG':'PG',
	'DYPG':'PG',
	'DOPG':'PG',
	'DGPG':'PG',
	'DEPG':'PG',
	'DNPG':'PG',
	'DLPS':'PS',
	'DMPS':'PS',
	'DPPS':'PS',
	'DSPS':'PS',
	'POPS':'PS',
	'PLPS':'PS',
	'SOPS':'PS',
	'SLPS':'PS',
	'DYPS':'PS',
	'YOPS':'PS',
	'DOPS':'PS',
	'DGPS':'PS',
	'DEPS':'PS',
	'DNPS':'PS',
	'DMPI':'PI',
	'DMPI13':'PI',
	'DMPI14':'PI',
	'DMPI15':'PI',
	'DMPI24':'PI',
	'DMPI25':'PI',
	'DMPI2A':'PI',
	'DMPI2B':'PI',
	'DMPI2C':'PI',
	'DMPI2D':'PI',
	'DMPI33':'PI',
	'DMPI34':'PI',
	'DMPI35':'PI',
	'PYPI':'PI',
	'POPI':'PI',
	'POPI13':'PI',
	'POPI14':'PI',
	'POPI15':'PI',
	'POPI24':'PI',
	'POPI25':'PI',
	'POPI2A':'PI',
	'POPI2B':'PI',
	'POPI2C':'PI',
	'POPI2D':'PI',
	'POPI33':'PI',
	'POPI34':'PI',
	'POPI35':'PI',
	'PLPI':'PI',
	'PLPI13':'PI',
	'PLPI14':'PI',
	'PLPI15':'PI',
	'PLPI24':'PI',
	'PLPI25':'PI',
	'PLPI2A':'PI',
	'PLPI2B':'PI',
	'PLPI2C':'PI',
	'PLPI2D':'PI',
	'PLPI33':'PI',
	'PLPI34':'PI',
	'PLPI35':'PI',
	'PNPI':'PI',
	'PNPI13':'PI',
	'PNPI14':'PI',
	'PNPI15':'PI',
	'PNPI24':'PI',
	'PNPI25':'PI',
	'PNPI2A':'PI',
	'PNPI2B':'PI',
	'PNPI2C':'PI',
	'PNPI2D':'PI',
	'PNPI33':'PI',
	'PNPI34':'PI',
	'PNPI35':'PI',
	'SAPI':'PI',
	'SAPI13':'PI',
	'SAPI14':'PI',
	'SAPI15':'PI',
	'SAPI24':'PI',
	'SAPI25':'PI',
	'SAPI2A':'PI',
	'SAPI2B':'PI',
	'SAPI2C':'PI',
	'SAPI2D':'PI',
	'SAPI33':'PI',
	'SAPI34':'PI',
	'SAPI35':'PI',
	'TMCL1':'CL',
	'TMCL2':'CL',
	'PMCL1':'CL',
	'PMCL2':'CL',
	'PVCL2':'CL',
	'TYCL1':'CL',
	'TYCL2':'CL',
	'TOCL1':'CL',
	'TOCL2':'CL',
	'LOACL1':'CL',
	'LOACL2':'CL',
	'LOCCL1':'CL',
	'LOCCL2':'CL',
	'TLCL1':'CL',
	'TLCL2':'CL',
	'LNCCL1':'CL',
	'LNCCL2':'CL',
	'LNACL1':'CL',
	'LNACL2':'CL',
	'LNDCL1':'CL',
	'LNDCL2':'CL',
	'LNBCL1':'CL',
	'LNBCL2':'CL',
	'SAPA':'PUFA',
	'SAPC':'PUFA',
	'SAPE':'PUFA',
	'SAPG':'PUFA',
	'SAPS':'PUFA',
	'SDPA':'PUFA',
	'SDPC':'PUFA',
	'SDPE':'PUFA',
	'SDPG':'PUFA',
	'SDPS':'PUFA',
	'DAPA':'PUFA',
	'DAPC':'PUFA',
	'DAPE':'PUFA',
	'DAPG':'PUFA',
	'DAPS':'PUFA',
	'PSM':'SM-N',
	'SSM':'SM-N',
	'ASM':'SM-N',
	'BSM':'SM-N',
	'23SM':'SM-N',
	'LSM':'SM-N',
	'OSM':'SM-N',
	'NSM':'SM-N',
	'CER160':'SM-P',
	'CER180':'SM-P',
	'CER181':'SM-P',
	'CER200':'SM-P',
	'CER220':'SM-P',
	'CER240':'SM-P',
	'CER241':'SM-P',
	'QMPE':'Bacterial Lipids',
	'PMPE':'Bacterial Lipids',
	'PMPG':'Bacterial Lipids',
	'PPPE':'Bacterial Lipids',
	'PVPE':'Bacterial Lipids',
	'PVPG':'Bacterial Lipids',
	'APPC':'Bacterial Lipids',
	'IPPC':'Bacterial Lipids',
	'PhPC':'Bacterial Lipids',
	'LAU':'Fatty Acids',
	'MYR':'Fatty Acids',
	'PAL':'Fatty Acids',
	'STE':'Fatty Acids',
	'ARA':'Fatty Acids',
	'BEH':'Fatty Acids',
	'TRI':'Fatty Acids',
	'LIGN':'Fatty Acids',
	'MYRO':'Fatty Acids',
	'PALO':'Fatty Acids',
	'HTA':'Fatty Acids',
	'OLE':'Fatty Acids',
	'LIN':'Fatty Acids',
	'ALIN':'Fatty Acids',
	'SDA':'Fatty Acids',
	'GLA':'Fatty Acids',
	'EICO':'Fatty Acids',
	'EDA':'Fatty Acids',
	'MEA':'Fatty Acids',
	'DGLA':'Fatty Acids',
	'ETE':'Fatty Acids',
	'ETA':'Fatty Acids',
	'EPA':'Fatty Acids',
	'ARAN':'Fatty Acids',
	'HPA':'Fatty Acids',
	'ERU':'Fatty Acids',
	'DDA':'Fatty Acids',
	'ADR':'Fatty Acids',
	'DPT':'Fatty Acids',
	'DPA':'Fatty Acids',
	'DHA':'Fatty Acids',
	'NER':'Fatty Acids',
	'TTA':'Fatty Acids',
	'TPT':'Fatty Acids',
	'TPA':'Fatty Acids',
	'THA':'Fatty Acids',
	'LAUP':'Fatty Acids',
	'MYRP':'Fatty Acids',
	'PALP':'Fatty Acids',
	'STEP':'Fatty Acids',
	'ARAP':'Fatty Acids',
	'BEHP':'Fatty Acids',
	'TRIP':'Fatty Acids',
	'LIGNP':'Fatty Acids',
	'MYROP':'Fatty Acids',
	'PALOP':'Fatty Acids',
	'HTAP':'Fatty Acids',
	'OLEP':'Fatty Acids',
	'LINP':'Fatty Acids',
	'ALINP':'Fatty Acids',
	'SDAP':'Fatty Acids',
	'GLAP':'Fatty Acids',
	'EICOP':'Fatty Acids',
	'EDAP':'Fatty Acids',
	'MEAP':'Fatty Acids',
	'DGLAP':'Fatty Acids',
	'ETEP':'Fatty Acids',
	'ETAP':'Fatty Acids',
	'EPAP':'Fatty Acids',
	'ARANP':'Fatty Acids',
	'HPAP':'Fatty Acids',
	'ERUP':'Fatty Acids',
	'DDAP':'Fatty Acids',
	'ADRP':'Fatty Acids',
	'DPTP':'Fatty Acids',
	'DPAP':'Fatty Acids',
	'DHAP':'Fatty Acids',
	'NERP':'Fatty Acids',
	'TTAP':'Fatty Acids',
	'TPTP':'Fatty Acids',
	'TPAP':'Fatty Acids',
	'THAP':'Fatty Acids',
	'SDS':'Detergents-P',
	'LMPG':'Detergents-P',
	'DHPC':'Detergents-P',
	'FOS10':'Detergents-P',
	'DPC':'Detergents-P',
	'TPC':'Detergents-P',
	'ADDG':'Detergents-2C1',
	'BDDG':'Detergents-2C1',
	'ADG':'Detergents-2C1',
	'BDG':'Detergents-2C1',
	'AOG':'Detergents-2C1',
	'BOG':'Detergents-2C1',
	'ADDM':'Detergents-2O1',
	'BDDM':'Detergents-2O1',
	'ADM':'Detergents-2O1',
	'BDM':'Detergents-2O1',
	'AOM':'Detergents-2O1',
	'BOM':'Detergents-2O1',
	'SB3-12':'Detergents-N',
	'SB3-14':'Detergents-N',
	'LDAO':'Detergents-N'
}

parent_atom = {
 	'Bacterial Lipids':'P',
	'CL':'C3',
	'Fatty Acids':'C1',
	'PA':'P',
	'PC':'P',
	'PE':'P',
	'PG':'P',
	'PI':'P',
	'PS':'P',
	'PUFA':'P',
	'Sterols':'O3',
	'SM-N':'N',
	'SM-NF':'NF',
	'Detergents-P':'P',
	'Detergents-N':'N',
	'Detergents-2C1':'2C1',
	'Detergents-2O1':'2O1'
}

lipidSurfArea = {
 	'Cholesterol':'40',
	'ERG':'55',
	'DLPA':'57.2',
	'DMPA':'63.5',
	'DPPA':'63',
	'DSPA':'61.8',
	'POPA':'60.1',
	'PLPA':'67.4',
	'SOPA':'58.8',
	'SLPA':'67.4',
	'DYPA':'67.4',
	'YOPA':'78',
	'DOPA':'64.8',
	'DGPA':'67.4',
	'DEPA':'67.4',
	'DNPA':'67.4',
	'DDPC':'64.1',
	'DCPC':'64.1',
	'DLPC':'64.1',
	'DMPC':'64.1',
	'DPPC':'63',
	'DSPC':'65.5',
	'POPC':'68.3',
	'PLPC':'65.2',
	'SOPC':'66.2',
	'SLPC':'65.2',
	'DYPC':'67',
	'YOPC':'67',
	'DOPC':'69.7',
	'DUPC':'70',
	'DGPC':'66.7',
	'DEPC':'66.5',
	'DNPC':'62.3',
	'DLPE':'60.8',
	'DMPE':'59.9',
	'DPPE':'59',
	'DSPE':'58.8',
	'PYPE':'58.8',
	'POPE':'58.8',
	'PLPE':'60.4',
	'SOPE':'56.9',
	'SLPE':'58.6',
	'DYPE':'62',
	'YOPE':'62',
	'OYPE':'61.7',
	'DOPE':'63.4',
	'DGPE':'61.4',
	'DEPE':'59',
	'DNPE':'57.5',
	'DLPG':'57.2',
	'DMPG':'60.6',
	'DPPG':'63',
	'DSPG':'63',
	'PYPG':'68.5',
	'POPG':'62',
	'PLPG':'67.4',
	'SOPG':'67.4',
	'SLPG':'67.4',
	'DYPG':'67.4',
	'DOPG':'67.4',
	'DGPG':'67.4',
	'DEPG':'67.4',
	'DNPG':'67.4',
	'DLPS':'58.1',
	'DMPS':'58.1',
	'DPPS':'60.8',
	'DSPS':'60',
	'POPS':'60.4',
	'PLPS':'61.7',
	'SOPS':'59.9',
	'SLPS':'60.7',
	'DYPS':'67.4',
	'YOPS':'65',
	'DOPS':'71.6',
	'DGPS':'65.1',
	'DEPS':'67.4',
	'DNPS':'67.4',
	'DMPI':'67.4',
	'DMPI13':'67.4',
	'DMPI14':'67.4',
	'DMPI15':'67.4',
	'DMPI24':'67.4',
	'DMPI25':'67.4',
	'DMPI2A':'67.4',
	'DMPI2B':'67.4',
	'DMPI2C':'67.4',
	'DMPI2D':'67.4',
	'DMPI33':'67.4',
	'DMPI34':'67.4',
	'DMPI35':'67.4',
	'PYPI':'67.4',
	'POPI':'67.4',
	'POPI13':'67.4',
	'POPI14':'67.4',
	'POPI15':'67.4',
	'POPI24':'67.4',
	'POPI25':'67.4',
	'POPI2A':'67.4',
	'POPI2B':'67.4',
	'POPI2C':'67.4',
	'POPI2D':'67.4',
	'POPI33':'67.4',
	'POPI34':'67.4',
	'POPI35':'67.4',
	'PLPI':'67.4',
	'PLPI13':'67.4',
	'PLPI14':'67.4',
	'PLPI15':'67.4',
	'PLPI24':'67.4',
	'PLPI25':'67.4',
	'PLPI2A':'67.4',
	'PLPI2B':'67.4',
	'PLPI2C':'67.4',
	'PLPI2D':'67.4',
	'PLPI33':'67.4',
	'PLPI34':'67.4',
	'PLPI35':'67.4',
	'PNPI':'67.4',
	'PNPI13':'67.4',
	'PNPI14':'67.4',
	'PNPI15':'67.4',
	'PNPI24':'67.4',
	'PNPI25':'67.4',
	'PNPI2A':'67.4',
	'PNPI2B':'67.4',
	'PNPI2C':'67.4',
	'PNPI2D':'67.4',
	'PNPI33':'67.4',
	'PNPI34':'67.4',
	'PNPI35':'67.4',
	'SAPI':'67.4',
	'SAPI13':'67.4',
	'SAPI14':'67.4',
	'SAPI15':'67.4',
	'SAPI24':'67.4',
	'SAPI25':'67.4',
	'SAPI2A':'67.4',
	'SAPI2B':'67.4',
	'SAPI2C':'67.4',
	'SAPI2D':'67.4',
	'SAPI33':'67.4',
	'SAPI34':'67.4',
	'SAPI35':'67.4',
	'TMCL1':'130',
	'TMCL2':'112.1',
	'PMCL1':'130',
	'PMCL2':'130',
	'PVCL2':'130',
	'TYCL1':'129',
	'TYCL2':'130.5',
	'TOCL1':'126.1',
	'TOCL2':'129.7',
	'LOACL1':'130',
	'LOACL2':'130',
	'LOCCL1':'130',
	'LOCCL2':'130',
	'TLCL1':'132.1',
	'TLCL2':'136.7',
	'LNCCL1':'130',
	'LNCCL2':'130',
	'LNACL1':'130',
	'LNACL2':'130',
	'LNDCL1':'130',
	'LNDCL2':'130',
	'LNBCL1':'130',
	'LNBCL2':'130',
	'SAPA':'66.2',
	'SAPC':'71.1',
	'SAPE':'64.7',
	'SAPG':'75.4',
	'SAPS':'65.9',
	'SDPA':'66.2',
	'SDPC':'67.2',
	'SDPE':'63.3',
	'SDPG':'67.4',
	'SDPS':'65.3',
	'DAPA':'72.8',
	'DAPC':'76.1',
	'DAPE':'70',
	'DAPG':'79.9',
	'DAPS':'71.7',
	'PSM':'55.4',
	'SSM':'55.3',
	'ASM':'55.4',
	'BSM':'55.4',
	'23SM':'55.4',
	'LSM':'55.4',
	'OSM':'55.4',
	'NSM':'55.4',
	'CER160':'55.3',
	'CER180':'55.3',
	'CER181':'55.3',
	'CER200':'55.3',
	'CER220':'55.3',
	'CER240':'55.3',
	'CER241':'55.3',
	'QMPE':'62.2',
	'PMPE':'60.5',
	'PMPG':'70.2',
	'PPPE':'63',
	'PVPE':'63',
	'PVPG':'62',
	'APPC':'61.1',
	'IPPC':'59.2',
	'PhPC':'79.5',
	'LAU':'30',
	'MYR':'30',
	'PAL':'30',
	'STE':'30',
	'ARA':'30',
	'BEH':'30',
	'TRI':'30',
	'LIGN':'30',
	'MYRO':'30',
	'PALO':'30',
	'HTA':'30',
	'OLE':'30',
	'LIN':'30',
	'ALIN':'30',
	'SDA':'30',
	'GLA':'30',
	'EICO':'30',
	'EDA':'30',
	'MEA':'30',
	'DGLA':'30',
	'ETE':'30',
	'ETA':'30',
	'EPA':'30',
	'ARAN':'30',
	'HPA':'30',
	'ERU':'30',
	'DDA':'30',
	'ADR':'30',
	'DPT':'30',
	'DPA':'30',
	'DHA':'30',
	'NER':'30',
	'TTA':'30',
	'TPT':'30',
	'TPA':'30',
	'THA':'30',
	'LAUP':'30',
	'MYRP':'30',
	'PALP':'30',
	'STEP':'30',
	'ARAP':'30',
	'BEHP':'30',
	'TRIP':'30',
	'LIGNP':'30',
	'MYROP':'30',
	'PALOP':'30',
	'HTAP':'30',
	'OLEP':'30',
	'LINP':'30',
	'ALINP':'30',
	'SDAP':'30',
	'GLAP':'30',
	'EICOP':'30',
	'EDAP':'30',
	'MEAP':'30',
	'DGLAP':'30',
	'ETEP':'30',
	'ETAP':'30',
	'EPAP':'30',
	'ARANP':'30',
	'HPAP':'30',
	'ERUP':'30',
	'DDAP':'30',
	'ADRP':'30',
	'DPTP':'30',
	'DPAP':'30',
	'DHAP':'30',
	'NERP':'30',
	'TTAP':'30',
	'TPTP':'30',
	'TPAP':'30',
	'THAP':'30',
	'SDS':'30',
	'LMPG':'30',
	'DHPC':'60',
	'FOS10':'30',
	'DPC':'30',
	'TPC':'30',
	'ADDG':'30',
	'BDDG':'30',
	'ADG':'30',
	'BDG':'30',
	'AOG':'30',
	'BOG':'30',
	'ADDM':'30',
	'BDDM':'30',
	'ADM':'30',
	'BDM':'30',
	'AOM':'30',
	'BOM':'30',
	'SB3-12':'30',
	'SB3-14':'30',
	'LDAO':'30'
}

membraneSize = {
 	'Cholesterol':'34.641',
	'ERG':'40.6202',
	'DLPA':'41.4246',
	'DMPA':'43.6463',
	'DPPA':'43.4741',
	'DSPA':'43.0581',
	'POPA':'42.4617',
	'PLPA':'44.9667',
	'SOPA':'42',
	'SLPA':'44.9667',
	'DYPA':'44.9667',
	'YOPA':'48.3735',
	'DOPA':'44.0908',
	'DGPA':'44.9667',
	'DEPA':'44.9667',
	'DNPA':'44.9667',
	'DDPC':'43.852',
	'DCPC':'43.852',
	'DLPC':'43.852',
	'DMPC':'43.852',
	'DPPC':'43.4741',
	'DSPC':'44.3283',
	'POPC':'45.2659',
	'PLPC':'44.2267',
	'SOPC':'44.5646',
	'SLPC':'44.2267',
	'DYPC':'44.833',
	'YOPC':'44.833',
	'DOPC':'45.7275',
	'DUPC':'45.8258',
	'DGPC':'44.7325',
	'DEPC':'44.6654',
	'DNPC':'43.2319',
	'DLPE':'42.7083',
	'DMPE':'42.391',
	'DPPE':'42.0714',
	'DSPE':'42',
	'PYPE':'42',
	'POPE':'42',
	'PLPE':'42.5676',
	'SOPE':'41.3159',
	'SLPE':'41.9285',
	'DYPE':'43.1277',
	'YOPE':'43.1277',
	'OYPE':'43.0232',
	'DOPE':'43.6119',
	'DGPE':'42.9185',
	'DEPE':'42.0714',
	'DNPE':'41.5331',
	'DLPG':'41.4246',
	'DMPG':'42.638',
	'DPPG':'43.4741',
	'DSPG':'43.4741',
	'PYPG':'45.3321',
	'POPG':'43.1277',
	'PLPG':'44.9667',
	'SOPG':'44.9667',
	'SLPG':'44.9667',
	'DYPG':'44.9667',
	'DOPG':'44.9667',
	'DGPG':'44.9667',
	'DEPG':'44.9667',
	'DNPG':'44.9667',
	'DLPS':'41.7493',
	'DMPS':'41.7493',
	'DPPS':'42.7083',
	'DSPS':'42.4264',
	'POPS':'42.5676',
	'PLPS':'43.0232',
	'SOPS':'42.391',
	'SLPS':'42.6732',
	'DYPS':'44.9667',
	'YOPS':'44.1588',
	'DOPS':'46.3465',
	'DGPS':'44.1928',
	'DEPS':'44.9667',
	'DNPS':'44.9667',
	'DMPI':'44.9667',
	'DMPI13':'44.9667',
	'DMPI14':'44.9667',
	'DMPI15':'44.9667',
	'DMPI24':'44.9667',
	'DMPI25':'44.9667',
	'DMPI2A':'44.9667',
	'DMPI2B':'44.9667',
	'DMPI2C':'44.9667',
	'DMPI2D':'44.9667',
	'DMPI33':'44.9667',
	'DMPI34':'44.9667',
	'DMPI35':'44.9667',
	'PYPI':'44.9667',
	'POPI':'44.9667',
	'POPI13':'44.9667',
	'POPI14':'44.9667',
	'POPI15':'44.9667',
	'POPI24':'44.9667',
	'POPI25':'44.9667',
	'POPI2A':'44.9667',
	'POPI2B':'44.9667',
	'POPI2C':'44.9667',
	'POPI2D':'44.9667',
	'POPI33':'44.9667',
	'POPI34':'44.9667',
	'POPI35':'44.9667',
	'PLPI':'44.9667',
	'PLPI13':'44.9667',
	'PLPI14':'44.9667',
	'PLPI15':'44.9667',
	'PLPI24':'44.9667',
	'PLPI25':'44.9667',
	'PLPI2A':'44.9667',
	'PLPI2B':'44.9667',
	'PLPI2C':'44.9667',
	'PLPI2D':'44.9667',
	'PLPI33':'44.9667',
	'PLPI34':'44.9667',
	'PLPI35':'44.9667',
	'PNPI':'44.9667',
	'PNPI13':'44.9667',
	'PNPI14':'44.9667',
	'PNPI15':'44.9667',
	'PNPI24':'44.9667',
	'PNPI25':'44.9667',
	'PNPI2A':'44.9667',
	'PNPI2B':'44.9667',
	'PNPI2C':'44.9667',
	'PNPI2D':'44.9667',
	'PNPI33':'44.9667',
	'PNPI34':'44.9667',
	'PNPI35':'44.9667',
	'SAPI':'44.9667',
	'SAPI13':'44.9667',
	'SAPI14':'44.9667',
	'SAPI15':'44.9667',
	'SAPI24':'44.9667',
	'SAPI25':'44.9667',
	'SAPI2A':'44.9667',
	'SAPI2B':'44.9667',
	'SAPI2C':'44.9667',
	'SAPI2D':'44.9667',
	'SAPI33':'44.9667',
	'SAPI34':'44.9667',
	'SAPI35':'44.9667',
	'TMCL1':'62.45',
	'TMCL2':'57.9914',
	'PMCL1':'62.45',
	'PMCL2':'62.45',
	'PVCL2':'62.45',
	'TYCL1':'62.2093',
	'TYCL2':'62.57',
	'TOCL1':'61.5061',
	'TOCL2':'62.3779',
	'LOACL1':'62.45',
	'LOACL2':'62.45',
	'LOCCL1':'62.45',
	'LOCCL2':'62.45',
	'TLCL1':'62.9524',
	'TLCL2':'64.0391',
	'LNCCL1':'62.45',
	'LNCCL2':'62.45',
	'LNACL1':'62.45',
	'LNACL2':'62.45',
	'LNDCL1':'62.45',
	'LNDCL2':'62.45',
	'LNBCL1':'62.45',
	'LNBCL2':'62.45',
	'SAPA':'44.5646',
	'SAPC':'46.1844',
	'SAPE':'44.0568',
	'SAPG':'47.5605',
	'SAPS':'44.4635',
	'SDPA':'44.5646',
	'SDPC':'44.8999',
	'SDPE':'43.5775',
	'SDPG':'44.9667',
	'SDPS':'44.2606',
	'DAPA':'46.7333',
	'DAPC':'47.7807',
	'DAPE':'45.8258',
	'DAPG':'48.9592',
	'DAPS':'46.3789',
	'PSM':'40.7676',
	'SSM':'40.7308',
	'ASM':'40.7676',
	'BSM':'40.7676',
	'23SM':'40.7676',
	'LSM':'40.7676',
	'OSM':'40.7676',
	'NSM':'40.7676',
	'CER160':'40.7308',
	'CER180':'40.7308',
	'CER181':'40.7308',
	'CER200':'40.7308',
	'CER220':'40.7308',
	'CER240':'40.7308',
	'CER241':'40.7308',
	'QMPE':'43.1972',
	'PMPE':'42.6028',
	'PMPG':'45.8912',
	'PPPE':'43.4741',
	'PVPE':'43.4741',
	'PVPG':'43.1277',
	'APPC':'42.8135',
	'IPPC':'42.1426',
	'PhPC':'48.8365',
	'LAU':'30',
	'MYR':'30',
	'PAL':'30',
	'STE':'30',
	'ARA':'30',
	'BEH':'30',
	'TRI':'30',
	'LIGN':'30',
	'MYRO':'30',
	'PALO':'30',
	'HTA':'30',
	'OLE':'30',
	'LIN':'30',
	'ALIN':'30',
	'SDA':'30',
	'GLA':'30',
	'EICO':'30',
	'EDA':'30',
	'MEA':'30',
	'DGLA':'30',
	'ETE':'30',
	'ETA':'30',
	'EPA':'30',
	'ARAN':'30',
	'HPA':'30',
	'ERU':'30',
	'DDA':'30',
	'ADR':'30',
	'DPT':'30',
	'DPA':'30',
	'DHA':'30',
	'NER':'30',
	'TTA':'30',
	'TPT':'30',
	'TPA':'30',
	'THA':'30',
	'LAUP':'30',
	'MYRP':'30',
	'PALP':'30',
	'STEP':'30',
	'ARAP':'30',
	'BEHP':'30',
	'TRIP':'30',
	'LIGNP':'30',
	'MYROP':'30',
	'PALOP':'30',
	'HTAP':'30',
	'OLEP':'30',
	'LINP':'30',
	'ALINP':'30',
	'SDAP':'30',
	'GLAP':'30',
	'EICOP':'30',
	'EDAP':'30',
	'MEAP':'30',
	'DGLAP':'30',
	'ETEP':'30',
	'ETAP':'30',
	'EPAP':'30',
	'ARANP':'30',
	'HPAP':'30',
	'ERUP':'30',
	'DDAP':'30',
	'ADRP':'30',
	'DPTP':'30',
	'DPAP':'30',
	'DHAP':'30',
	'NERP':'30',
	'TTAP':'30',
	'TPTP':'30',
	'TPAP':'30',
	'THAP':'30',
	'SDS':'30',
	'LMPG':'30',
	'DHPC':'42.4264',
	'FOS10':'30',
	'DPC':'30',
	'TPC':'30',
	'ADDG':'30',
	'BDDG':'30',
	'ADG':'30',
	'BDG':'30',
	'AOG':'30',
	'BOG':'30',
	'ADDM':'30',
	'BDDM':'30',
	'ADM':'30',
	'BDM':'30',
	'AOM':'30',
	'BOM':'30',
	'SB3-12':'30',
	'SB3-14':'30',
	'LDAO':'30'
}
