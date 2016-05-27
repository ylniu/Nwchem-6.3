      subroutine ebt_filename(file,dirfile)
      implicit none
#include "inp.fh"
#include "global.fh"
      character*(*) file,dirfile
c
      character*255 dirss
      integer lenn0,lenn1,i
c

      call util_directory_name(dirss, .false., ga_nodeid())
      dirfile(1:inp_strlen(dirss))=dirss(1:inp_strlen(dirss))
      dirfile(inp_strlen(dirss)+1:inp_strlen(dirss)+1)='/'
      lenn0=inp_strlen(dirss)+2
      lenn1=inp_strlen(dirss)+inp_strlen(file)+1
      dirfile(lenn0:lenn1)=file(1:inp_strlen(file))
      do i=lenn1+1,255
      dirfile(i:i)=' '
      enddo
      
      return
      end
