#!/usr/bin/env python3
import re, csv, sys, getopt, getpass, mysql.connector, argparse

def create_and_parse_argument_options(argument_list):
	parser = argparse.ArgumentParser(description='Update structural folds tables from ECOD and SCOPe, given a PDB ID')
	parser.add_argument('pdb_id', help='PDB identifier to query', type=str)
	parser.add_argument('user_name', help='Username for connecting to DESIRE', type=str)
	parser.add_argument('-host','--db_host', help='Defines database host (default: 130.207.36.75)', type=str, default='130.207.36.75')
	parser.add_argument('-schema','--db_schema', help='Defines schema to use (default: SEREB)', type=str, default='SEREB')
	parser.add_argument('-dl','--download_most_recent_fold_definitions', help='Update latest fold definitions.', default=False, action="store_true")
	commandline_args = parser.parse_args(argument_list)
	return commandline_args

def download_latest_fold_defs(url, file_name):
	'''Download latest structural phylogeny file'''
	import urllib.request
	try:
		urllib.request.urlretrieve(url, file_name)
	except:
		raise ConnectionError("Failed downloading "+ url)
	return True

def initiate_connection(uname, host, database):
	pw = getpass.getpass("Password: ")
	cnx = mysql.connector.connect(user=uname, password=pw, host=host, database=database)
	return cnx

def parse_ecod_definitions(pdbid):
	ecod_defs = list()
	with open("ecod.latest.domains.txt") as file:
		for line in file:
			if re.match(r'^#', line):
				continue
			if pdbid == line.split('\t')[4]:
				ecod_defs.append([line.split('\t')[3],line.split('\t')[5],line.split('\t')[7],line.split('\t')[8],
				line.split('\t')[9],line.split('\t')[10],line.split('\t')[11],line.split('\t')[12],line.split('\t')[13]])
	if len(ecod_defs) == 0:
		raise ValueError("Entered pdb_id: "+ pdbid+" was not found in the latest ecod file!")
	return ecod_defs

def upload_struc_fold(level, name, class_sys, parent, external_id, cursor, cnx):
	if level != 'Architecture':
		query = "INSERT INTO `SEREB`.`Structural_Folds`(`Level`,`Name`,`classification_system`,`parent`,`external_id`) VALUES('"+level+"','"+name+"','"+class_sys+"','"+parent+"','"+external_id+"')"
	else:
		query = "INSERT INTO `SEREB`.`Structural_Folds`(`Level`,`Name`,`classification_system`,`external_id`) VALUES('"+level+"','"+name+"','"+class_sys+"','"+external_id+"')"
	cursor.execute(query)
	lastrow_id = str(cursor.lastrowid)
	cnx.commit()
	return lastrow_id

def upload_strucfold_chains(cursor, strucfoldid, chainid):
	cursor.execute("SELECT * FROM SEREB.StrucFold_Chains WHERE\
					strucfold_id = '"+strucfoldid+"' AND\
					chain_id = '"+chainid+"'")
	result = cursor.fetchall()
	if len(result) == 0:
		query = "INSERT INTO `SEREB`.`StrucFold_Chains`(`strucfold_id`, `chain_id`) VALUES ('"+strucfoldid+"','"+chainid+"')"
		cursor.execute(query)
	return True

def upload_strucfold_resis(cursor, strucfoldid, resi_ids):
	for strucid, resid in map(lambda e: (e, strucfoldid), resi_ids):
		cursor.execute("SELECT * FROM SEREB.StrucFold_Residues WHERE\
						strucfold_id = '"+str(strucid)+"' AND\
						residue_id = '"+str(resid)+"'")
		result = cursor.fetchall()
		if len(result) == 0:
			query = "INSERT INTO `SEREB`.`StrucFold_Residues`(`residue_id`, `strucfold_id`) VALUES ('"+str(strucid)+"','"+str(resid)+"')"
			cursor.execute(query)
	return True

def get_resis(cursor, pol_id, resi_ranges):
	between_statement = ""
	if len(resi_ranges.split(",")) == 1:
		resirange1 = resi_ranges.split(":")[1].split("-")[0]
		resirange2 = resi_ranges.split(":")[1].split("-")[1]
		between_statement = "SEREB.Residues.resNum BETWEEN "+resirange1+" AND "+resirange2
	if len(resi_ranges.split(",")) > 1:
		list_betweens = list()
		for one_range in resi_ranges.split(","):
			resirange1 = one_range.split(":")[1].split("-")[0]
			resirange2 = one_range.split(":")[1].split("-")[1]
			list_betweens.append("SEREB.Residues.resNum BETWEEN "+resirange1+" AND "+resirange2)
		between_statement = " OR ".join(list_betweens)
	
	cursor.execute("SELECT resi_id FROM SEREB.Residues WHERE\
					SEREB.Residues.PolData_id = "+pol_id+" AND ("+between_statement+")")
	result = cursor.fetchall()
	return [item for sublist in result for item in sublist]

def check_polymerid_from_chainid(cursor, pdbid, chainid):
	cursor.execute("SELECT polymer_id, ChainList_id FROM SEREB.ChainList\
					INNER JOIN SEREB.ThreeDStructures ON SEREB.ChainList.3D_structure_id = SEREB.ThreeDStructures.3D_structure_id\
					WHERE SEREB.ThreeDStructures.StructureName = '"+pdbid.upper()+"' AND\
					SEREB.ChainList.ChainName = '"+chainid+"'")
	result = cursor.fetchall()
	if len(result) == 0:
		return None, None
	if len(result) == 1:
		return result[0][0], result[0][1]
	if len(result) > 1:
		raise ValueError("Mistake in database for query with pdb "+pdbid+" and chain "+chainid+"!")

def check_then_upload_struc_fold(ecod_definitions, cursor, cnx, pdbid):
	levels = ['Architecture', 'X', 'H', 'T', 'F']
	for entry_for_upload in ecod_definitions:
		if re.match(r'UNCLASSIFIED', entry_for_upload[8]):
			continue
		parent = "NULL"
		last_statement = "IS "+parent
		external_ecod_id = ['0'] + entry_for_upload[0].split('.')
		###------Separate this into function
		for fold_name,fold_level,fold_id in zip(entry_for_upload[-5:], levels, external_ecod_id):
			fold_name = fold_name.replace("\"","")
			if fold_level != 'Architecture':
				last_statement = "= '"+parent+"'"
			cursor.execute("SELECT SEREB.Structural_Folds.struc_fold_id FROM SEREB.Structural_Folds WHERE\
							SEREB.Structural_Folds.Level = '"+fold_level+"' AND\
							SEREB.Structural_Folds.Name = '"+fold_name+"' AND\
							SEREB.Structural_Folds.external_id = '"+fold_id+"' AND\
							SEREB.Structural_Folds.parent "+last_statement)
			result = cursor.fetchall()
			if len(result) == 0:
				parent = upload_struc_fold(fold_level, fold_name, 'ECOD', parent, fold_id, cursor, cnx)
			if len(result) == 1:
				parent = str(result[0][0])
		###------Up to here
		pol_id, chain_id = check_polymerid_from_chainid(cursor, pdbid, entry_for_upload[1])
		if chain_id is not None:
			upload_strucfold_chains(cursor, str(parent), str(chain_id))
		if pol_id is not None:
			resi_ids_for_upload = get_resis(cursor, str(pol_id), entry_for_upload[2])
			if len(resi_ids_for_upload) > 0:
				upload_strucfold_resis(cursor, str(parent), resi_ids_for_upload)

def main(commandline_arguments):
	comm_args = create_and_parse_argument_options(commandline_arguments)
	pdbid = str(comm_args.pdb_id).lower()
	if comm_args.download_most_recent_fold_definitions:
		download_latest_fold_defs('http://prodata.swmed.edu/ecod/distributions/ecod.latest.domains.txt', 'ecod.latest.domains.txt')
	ecod_defs = parse_ecod_definitions(pdbid)

	cnx = initiate_connection(comm_args.user_name, comm_args.db_host, comm_args.db_schema)
	cursor = cnx.cursor()

	check_then_upload_struc_fold(ecod_defs, cursor, cnx, pdbid)
	
	cnx.commit()
	cursor.close()
	cnx.close()

if __name__ == '__main__':
	sys.exit(main(sys.argv[1:]))