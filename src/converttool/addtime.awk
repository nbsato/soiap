/^!#TIMEMOD/{
  print "      use m_showtime"
  next
}
/^!#TIMEDEF/{
  print "      type(t_showtime)::  timer"
  next
}
/^!#TIME0/{
  print "      call timer%start()"
  next
}
/^!#TIME1/{
  print "      call timer%stop()"
  next
}
/^!#TIMESHOW/{
  s=$2
  print "      call timer%show("s")"
  next
}
{ print}
