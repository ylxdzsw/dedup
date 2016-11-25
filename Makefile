PYEXECUTABLE = python

.PHONY: c_ext clean

c_ext: compile.py
	$(PYEXECUTABLE) $<

clean:
	rm -f *.so *.o
