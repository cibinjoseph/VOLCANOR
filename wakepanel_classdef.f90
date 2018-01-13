module wakepanel_classdef
  use vr_classdef
  implicit none
  type wakepanel_class
    type(vr_class) :: vr
    integer :: tag                   ! for identifying panel to be wing or wake
  end type wakepanel_class

contains

  ! VR coordinates
  ! o---------> Y along span
  ! |
  ! |   1-----------4
  ! |   |     4     |
  ! |   |           |
  ! |   |1         3|
  ! |   |           |
  ! |   |     2     |
  ! |   2-----------3
  ! |
  ! V X along chord

end module wakepanel_classdef
