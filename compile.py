
import sys
import os
import inspect

cur_dir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
print cur_dir, '<<'

print 'compiling Infomap'
os.chdir('./Sources/Infomap-0.11.5/')
os.system('make clean')
os.system('make')
#print 'leaving directory: ',
#os.system('pwd')
os.chdir('../..')
#print 'moved back to directory: ',
#os.system('pwd')
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

print '\n\n\n============================================== '
print 'Running code with no arguments: '
os.system('./bin/topicmap')

print '\n\n\n============================================== '
print 'The code has been installed in ', cur_dir
print 'Please do not rename/move the folder "bin" and its content from', cur_dir
print 'Since topicmap calls ', infomap_path, ', renaming the executables breaks the code.'
print 'If you are not happy with the current location,'
print 'please move the entire folder "topicmapping" to the desired location'
print 'and run this script again.'













