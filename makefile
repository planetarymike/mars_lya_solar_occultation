IDIR=-I./
CC=g++
LIBS=$(shell pkg-config --cflags --libs gsl)
SRCFNSLOCFLG=-D 'SRCFNSLOC="./"'

solar_LOS_profiles:
	$(CC) compute_solar_LOS_profiles.cpp $(IDIR) $(LIBS) $(SRCFNSLOCFLG) -g -O3 -o compute_solar_LOS_profiles.x
