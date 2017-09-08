module pymon
contains
FUNCTION getnames() result(out_names)
    use dsmacc_parameters, only : nspec
    use dsmacc_monitor, only : spc_names
    !integer, parameter : nspec = 611
    character, dimension(nspec,24) :: out_names
    integer i,j
    DO i=1,NSPEC
      do j=1,24
          out_names(i,j) = spc_names(i)(j:j)
      end do
    END DO
end FUNCTION
end module pymon