#! /usr/bin/env python
from ete3 import NCBITaxa
import sys,os,argparse,shutil
def main(argv):



    parser = argparse.ArgumentParser(description="Initialize ncbi taxonomy DB without overwriting anything currently existing")
    parser.add_argument("--taxdump", type=str, help="Path to the taxdump.tar.gz file")
    parser.add_argument("--dbpath", type=str, help="Path to the directory where the database should be saved")

    args = parser.parse_args()

    taxdump_path = args.taxdump
    dbpath = args.dbpath
    homedir=os.path.expanduser('~')
    if (os.path.exists(homedir + '/.etetoolkit/') and os.path.exists(homedir + '/.etetoolkit/taxa.sqlite')):
       #copy it
       print("Copying old taxa.sqlite files to tmp files")
       os.makedirs(homedir + '/.etetoolkit/tmp')
       shutil.copy(homedir + '/.etetoolkit/taxa.sqlite', homedir + '/.etetoolkit/tmp/taxa.sqlite')
       shutil.copy(homedir + '/.etetoolkit/taxa.sqlite.traverse.pkl', homedir + '/.etetoolkit/tmp/taxa.sqlite.traverse.pkl')

    #init

    print("Creating database. If you have never initialized NCBI taxonomy database before, this will take longer - first, ete3 will download the current NCBI database and then your own taxonomy database will be used to overwrite it.")
    ncbi = NCBITaxa()
    
    ncbi.update_taxonomy_database(taxdump_file=taxdump_path)

    #move to new dir
    print("Moving database to " + dbpath)
    shutil.copy(homedir + '/.etetoolkit/taxa.sqlite', dbpath+"/taxa.sqlite")
    shutil.copy(homedir + '/.etetoolkit/taxa.sqlite.traverse.pkl', dbpath + "/taxa.sqlite.traverse.pkl")
    os.remove(homedir + '/.etetoolkit/taxa.sqlite')
    os.remove(homedir + '/.etetoolkit/taxa.sqlite.traverse.pkl')

    #if tmpdir exists move back
    if os.path.exists(homedir + '/.etetoolkit/tmp/taxa.sqlite'):
        print("Restoring tmp files")
        os.rename(homedir + '/.etetoolkit/tmp/taxa.sqlite',homedir + '/.etetoolkit/taxa.sqlite')
        os.rename(homedir + '/.etetoolkit/tmp/taxa.sqlite.traverse.pkl', homedir + '/.etetoolkit/taxa.sqlite.traverse.pkl')
        os.rmdir(homedir + '/.etetoolkit/tmp/')
    
if __name__ == "__main__":
  main(sys.argv)
