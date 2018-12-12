c-----
c       handler to catch floating point errros for SUN
c       do not use in production work, since it uses a
c       lot of cycles
c-----
c       do a f77 -g compile, and
c       use dbx, catch FPE, run
c       the dbx  will stop at floating point errors
c-----
      subroutine mchdep()
      external handler
      call ieee_handler('set','common',handler)
      return
      end

      subroutine handler(i1,i2,i3,i4)
      return
      end
