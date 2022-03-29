EIGEN=$(EIGEN_DIR)
CONEX=/home/frank/conexnew/conex

main.o: main.cc
	g++ -c -I$(EIGEN) -I$(CONEX) -g --std=c++17  main.cc

a.out:  main.o
	g++ --std=c++17  main.o  


