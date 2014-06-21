
import sys
import os
import inspect

cur_dir= os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
print cur_dir, '<<'

print 'compiling Infomap'
os.chdir('./Sources/Infomap-0.15.7/')
os.system('make clean')
os.system('make')
if os.path.isfile('Infomap') == False:
    print 'Infomap was not compiled. Trying again...'
    os.system('make clean')
    os.system('make noomp')
if os.path.isfile('Infomap') == False:
    print 'Infomap did not compile.'
    print 'Please contact me: arg.lanci@gmail.com'
    exit()

#print 'leaving directory: ',
#os.system('pwd')
os.chdir('../..')
#print 'moved back to directory: ',
#os.system('pwd')
os.system('rm -r bin/')
os.system('mkdir bin')
os.system('mv ./Sources/Infomap-0.15.7/Infomap bin/Infomap')


print 'compiling alpha_optimization'
os.system('g++ -O3 -funroll-loops -Wall -o bin/opt_alpha ./Sources/TopicMapping/optimize_alpha.cpp')

print 'compiling pajek formatter'
os.system('g++ -O3 -funroll-loops -Wall -o bin/edges2pajek ./Sources/SingleSliceInfomap/pajek_format.cpp')

print 'compiling partition converter'
os.system('g++ -O3 -funroll-loops -Wall -o bin/tree2part ./Sources/SingleSliceInfomap/get_partition.cpp')

print 'compiling show topic info'
os.system('g++ -O3 -funroll-loops -Wall -o bin/show_topic_info ./Sources/TopicMapping/show_topic_info.cpp')



infomap_path= cur_dir+'/bin/Infomap'

original_file='./Sources/TopicMapping/docmap.cpp'

replaced_file=open('./Sources/TopicMapping/docmap_tmp.cpp', 'w')

# asserting that there no empty spaces (Infomap's call would not work)
if len(infomap_path.split())!=1:
    print '\n\nERROR with directory name', cur_dir
    print 'PLEASE REMOVE EMPTY SPACES from directory name and compile again'
    exit(-1)

for l in open(original_file):
    replaced_file.write(l.replace('"___INSERT_INFOMAP_PATH_HERE___"', '"'+infomap_path+'"'))
replaced_file.close()


print 'compiling topic mapping'
os.system('g++ -O3 -funroll-loops -Wall -o bin/topicmap ./Sources/TopicMapping/docmap_tmp.cpp')

print '\n\n\n============================================== '
print 'Running code with no arguments: '
os.system('./bin/topicmap')

print '\n\n\n============================================== '
print 'This program has been installed in:'
print cur_dir
print 'Please do not rename/move the folder "bin" and its content from there.'
print 'Since topicmap calls:'
print infomap_path
print 'renaming the executables breaks the code.'
print 'If you are not happy with the current location,'
print 'please move the entire current folder to the desired location',
print 'and run this script again.'













