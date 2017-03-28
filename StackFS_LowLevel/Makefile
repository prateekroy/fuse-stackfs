all: StackFS_ll

StackFS_ll:
	g++ -Wall -Werror -D_FILE_OFFSET_BITS=64 StackFS_LowLevel.c `pkg-config fuse3 --cflags --libs glib-2.0` -o StackFS_ll
	cd leon/build/ make -j 4
clean:
	rm -rf StackFS_ll
