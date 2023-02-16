# $(addsuffix <prefix>, <name1 name2 ...>)
# add prefix as suffix in all names
# $(basename <names..>)
# fetch suffix in all names
# $(notdir <names..>)
# remove all documents path
# the 3 commands can generate obj names.

# $(patsubst <pattern>, <replacement>, <text>)
# if names in text match <pattern>, then replace it with <replacement>
# $(wildcard <pattern>)
# return all names match <pattern>
EXE = BRayTracer
CXX = g++
DESTDIR = ./

# project cpps
SOURCES += $(wildcard *.cpp)

OBJECTS = $(addsuffix .o, $(basename $(notdir $(SOURCES))))

LIBS     =
CXXFLAGS = -ggdb
LDFLAGS = 

##---------------------------------------------------------------------
## BUILD RULES
##---------------------------------------------------------------------

%.o:%.cpp
		$(CXX) $(CXXFLAGS) -c -o $@ $<

all: $(DESTDIR)$(EXE)
		@echo Build complete for Linux

$(EXE): $(OBJECTS)
		$(CXX) -o $@ $^ $(CXXFLAGS) $(LIBS)

clean:
		rm -f $(EXE) $(OBJECTS)

