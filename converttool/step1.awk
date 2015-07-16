/^ *subroutine input/{
print "subroutine show_argument()\n\
  implicit none\n\
  write(6,*)\n\
  write(6,*)'argument:'\n\
  write(6,*)'thisprogram inputfilename'\n\
  write(6,*)\n\
end subroutine show_argument"
}
/! unit cell/{
print "  character(120):: inputfilename\n\
  integer:: narg\n\
 narg=command_argument_count()\n\
  if (narg/=1) then \n\
     call show_argument()\n\
     stop 100\n\
  endif\n\
  call get_command_argument(1,inputfilename)\n\
  write(6,*)'input=',trim(inputfilename)"
}
/\"input.dat\"/{
 gsub("\"input\.dat\"","inputfilename")
}
{ print }
