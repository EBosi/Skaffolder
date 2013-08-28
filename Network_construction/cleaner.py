import os,sys,shutil

def check(dir_):
	if dir_[-1]=='/': return dir_[:-1]
	else: return dir_

inp=check(sys.argv[1])
backup='/'.join(inp.split('/')[:-1])+'/BACKUP/'
print backup
for file_ in os.listdir(inp):
	os.remove(inp+'/'+file_)
for file_ in os.listdir(backup):
	shutil.copy2(backup+'/'+file_,inp)
