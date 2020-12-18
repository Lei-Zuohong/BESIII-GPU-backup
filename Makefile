#========================================================
#       Multi-Architecture makefile...
#========================================================

ifeq (,$(filter _%,$(notdir $(CURDIR))))
include $(GPUPWA)/target.mk
else
ifeq (_common, $(notdir $(CURDIR)))

VPATH = $(SRCDIR)

.DEFAULT: ; @:


else

VPATH = $(SRCDIR):$(SRCDIR)/_common

.SUFFIXES : .o .c .C .h .cl .cpp .cxx


include $(GPUPWA)/paths.mk

include $(GPUPWA)/flags.mk


all: depend binfiles pipipi 

binfiles:
		@if [[ -f $(SRCDIR)/binfiles ]]; \
		then \
		  echo "binfile directory exists"; \
		else \
		  echo "Linking binfile directory"; \
		  ln -s $(GPUPWA)/GPUPWA/_common/binfiles $(SRCDIR)/binfiles;\
		fi


include $(GPUPWA)/depends.mk


PIPIPI_OBJS= $(GPUPWALIB)

pipipi: pipipi.o $(GPUPWALIB)
		$(CC) pipipi.o $(LDFLAGS) $(PIPIPI_OBJS) \
		-o pipipi 


$(GPUPWALIB): $(GPUPWADIR)/*.h $(GPUPWADIR)/*.cpp $(GPUPWADIR)/*.cl
		make -C $(GPUPWADIR) lib

include $(GPUPWA)/commands.mk

endif
endif
