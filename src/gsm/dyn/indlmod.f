! input forcing parameter module
      module indlmod

      implicit none

      contains

      subroutine get_indl(jba, n, l, indlsev)
        integer, intent(in)  :: jba, n, l
        integer, intent(out) :: indl

        indl = jba + (n-l)/2 + 1
      end subroutine get_indl

      end module indlmod
