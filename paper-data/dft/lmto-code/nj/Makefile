srcpath := $(dir $(firstword $(MAKEFILE_LIST)))
srcpath := $(srcpath:/=)
srcpath ?= .


# $(info $(srcpath))
# $(info $(MAKEFILE_LIST))

MAKEFLAGS += -r
MAKEFLAGS += -R
goals := $(filter-out nj clean browse commands deps graph query targets compdb recompact toollist stats explain keeprsp debuglist help,$(MAKECMDGOALS))
# $(info $(MAKECMDGOALS))
# $(info $(goal))

.PHONY: nj

build.py = $(srcpath)/build.py

$(goals) nj: | build.ninja
	ninja $(goals)

build.ninja: $(build.py)
	python $(build.py) $@

clean browse commands deps graph query targets compdb recompact toollist:
	ninja -t $(subst toollist,list,$@)

stats explain keeprsp debuglist:
	ninja -d $(subst debuglist,list,$@)

help:
	-@$(shell ninja -h) # to avoid having "make: [help] Error 1 (ignored)"

.DEFAULT_GOAL := nj
.SUFFIXES:

