interface Profiler_start
   subroutine Profiler_startName(name)
     character (len=*), intent(in) :: name
   end subroutine Profiler_startName

   subroutine Profiler_startId(id)
     integer, intent(in) :: id
  end subroutine Profiler_startId
end interface

interface Profiler_stop
   subroutine Profiler_stopName(name)
     character (len=*), intent(in) :: name
   end subroutine Profiler_stopName

   subroutine Profiler_stopId(id)
     integer, intent(in) :: id
   end subroutine Profiler_stopId
end interface


