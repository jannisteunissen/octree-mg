LIB := libfishpack.a
F77 := gfortran -O2 -std=legacy -fdefault-real-8 -fdefault-double-8
CPP := gfortran -cpp

SRC := blktri.f cblktri.f cmgnbn.f comf.f fftpack.f genbun.f gnbnaux.f hstcrt.f	\
    hstcsp.f hstcyl.f hstplr.f hstssp.f hw3crt.f hwscrt.f hwscsp.f hwscyl.f	\
    hwsplr.f hwsssp.f pois3d.f poistg.f sepaux.f sepeli.f sepx4.f

OBJS := $(subst .f,.o,$(SRC))

$(LIB) : $(OBJS)
	$(RM) $@
	$(AR) rcs $@ $^

%.o : %.f
	$(F77) -c $< -o $@

clean:
	rm -f $(LIB) $(OBJS) 
