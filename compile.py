
import sys
import os
import inspect

cur_dir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
print cur_dir, '<<'

os.system('pwd')
print 'compiling Infomap'
os.chdir('./Sources/Infomap-0.11.5/')
os.system('make clean')
os.system('make')
os.system('pwd')
os.chdir('../..')
os.system('pwd')
os.system('rm -r bin/')
os.system('mkdir bin')
os.system('mv ./Sources/Infomap-0.11.5/Infomap bin/Infomap')

infomap_path= cur_dir+'/bin/Infomap'

original_file='./Sources/TopicMapping/docmap.cpp'

replaced_file=open('./Sources/TopicMapping/docmap_tmp.cpp', 'w')

for l in open(original_file):
    replaced_file.write(l.replace('"___INSERT_INFOMAP_PATH_HERE___"', '"'+infomap_path+'"'))
replaced_file.close()


print 'compiling topic mapping'
os.system('g++ -O3 -funroll-loops -Wall -o bin/topicmap ./Sources/TopicMapping/docmap_tmp.cpp')










