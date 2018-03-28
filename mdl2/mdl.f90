! This file is part of PKDGRAV3 (http://www.pkdgrav.org/).
! Copyright (c) 2001-2018 Joachim Stadel & Douglas Potter
!
! PKDGRAV3 is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! PKDGRAV3 is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with PKDGRAV3.  If not, see <http://www.gnu.org/licenses/>.

module mdl
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int32_t, c_int16_t
  type,bind(c) :: mdl2_t
    integer(c_int32_t)::ncpu     ! Global number of threads (total)
    integer(c_int32_t)::myid     ! Global index of this thread
    integer(c_int32_t)::nrank;   ! Number of global processes (e.g., MPI ranks)
    integer(c_int32_t)::irank;   ! Index of this process (MPI rank)
    integer(c_int16_t)::ncore;   ! Number of threads in this process
    integer(c_int16_t)::icore;   ! Local core id
  end type mdl2_t

  interface
    subroutine mdl_service(ctx,In,nIn,Out,nOut) BIND(C)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
      implicit none
      type(*)::ctx
      integer(c_int),VALUE::nIn ! Size of the input structure
      integer(c_int)::nOut    ! Size of the result (set)
      type(*)::In    ! Input type/structure
      type(*)::Out   ! Output type/structure
    end subroutine mdl_service

    subroutine mdl_launch(argc,argv,master,worker) BIND(C,NAME="mdlLaunch")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR, C_FUNPTR
      INTEGER(KIND=C_INT), VALUE :: argc
      TYPE(C_PTR), VALUE :: argv
      type(c_funptr), intent(in), value :: master
      type(c_funptr), intent(in), value :: worker
    end subroutine mdl_launch

    subroutine mdl_add_service(mdl,sid,p1,service,nin,nout) BIND(C,NAME="mdlAddService")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FUNPTR
      import :: mdl2_t
      type(mdl2_t)::mdl
      INTEGER(C_INT), VALUE :: sid
      type(*) :: p1
      type(c_funptr), intent(in), value :: service
!      procedure(mdl_service),pointer,intent(in)::service
      integer(c_int), VALUE :: nin, nout
    end subroutine mdl_add_service

    subroutine mdl_commit_services(mdl) BIND(C,NAME="mdlCommitServices")
      import :: mdl2_t
      type(mdl2_t)::mdl
    end subroutine mdl_commit_services

    subroutine mdl_handler(mdl) BIND(C,NAME="mdlHandler")
      import :: mdl2_t
      type(mdl2_t)::mdl
    end subroutine mdl_handler

    integer function mdl_req_service(mdl,who,what,msg,n) BIND(C,NAME="mdlReqService")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
      import :: mdl2_t
      type(mdl2_t)::mdl
      type(*) :: msg
      integer(c_int), VALUE :: who, what, n
    end function mdl_req_service

    subroutine mdl_get_reply(mdl,rid,msg,n) BIND(C,NAME="mdlGetReply")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
      import :: mdl2_t
      type(mdl2_t)::mdl
      integer(c_int), VALUE :: rid
      type(*) :: msg
      integer(c_int) :: n
    end subroutine mdl_get_reply

    integer function mdl_proc_to_thread(mdl,iProc) BIND(C,NAME="mdlBaseProcToThread")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
      import :: mdl2_t
      type(mdl2_t)::mdl
      integer(c_int), VALUE :: iProc
    end function mdl_proc_to_thread

    integer function mdl_thread_to_proc(mdl,iThread) BIND(C,NAME="mdlBaseThreadToProc")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
      import :: mdl2_t
      type(mdl2_t)::mdl
      integer(c_int), VALUE :: iThread
    end function mdl_thread_to_proc

  end interface

  contains
    integer function mdl_threads(mdl)
      type(mdl2_t)::mdl
      mdl_threads = mdl%ncpu
    end function mdl_threads

    integer function mdl_self(mdl)
      type(mdl2_t)::mdl
      mdl_self = mdl%myid
    end function mdl_self

    integer function mdl_core(mdl)
      type(mdl2_t)::mdl
      mdl_core = mdl%icore
    end function mdl_core

    integer function mdl_cores(mdl)
      type(mdl2_t)::mdl
      mdl_cores = mdl%ncore
    end function mdl_cores

    integer function mdl_proc(mdl)
      type(mdl2_t)::mdl
      mdl_proc = mdl%irank
    end function mdl_proc

    integer function mdl_procs(mdl)
      type(mdl2_t)::mdl
      mdl_procs = mdl%nrank
    end function mdl_procs

end module mdl
