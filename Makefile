CC=g++
LOP=-o
LOPT=-O3 -funroll-loops -Wall

MAIN=./Sources/TopicMapping/docmap.cpp
TAG=topicmap


$(MAIN).o :
	$(CC) $(LOPT) $(LOP) $(TAG) $(MAIN)


