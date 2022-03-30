EIGEN=$(EIGEN_DIR)
CONEX=/home/frank/conexnew/conex

examples.o: example_graphs.cc
	g++ -c -I$(CONEX) -g --std=c++17  example_graphs.cc

main.o: main.cc
	g++ -c -I$(EIGEN) -I$(CONEX) -g --std=c++17  main.cc

a.out:  main.o example_graphs.o
	g++ --std=c++17  main.o  example_graphs.o


