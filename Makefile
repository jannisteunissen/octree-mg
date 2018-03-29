.PHONY:	all clean

all:
	$(MAKE) -C lib_2d
	$(MAKE) -C lib_3d

clean:
	$(MAKE) -C lib_2d clean
	$(MAKE) -C lib_3d clean
