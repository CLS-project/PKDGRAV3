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

module mdl_module
  USE, INTRINSIC :: ISO_C_BINDING, ONLY: c_int32_t, c_int16_t
  type,bind(c) :: mdl_t
    integer(c_int32_t)::ncpu     ! Global number of threads (total)
    integer(c_int32_t)::myid     ! Global index of this thread
    integer(c_int32_t)::nrank;   ! Number of global processes (e.g., MPI ranks)
    integer(c_int32_t)::irank;   ! Index of this process (MPI rank)
    integer(c_int16_t)::ncore;   ! Number of threads in this process
    integer(c_int16_t)::icore;   ! Local core id
  end type mdl_t

  interface mdl_get_reply
    module procedure mdl_get_reply_array, mdl_get_reply_scalar
  end interface mdl_get_reply

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

    subroutine mdl_launch(argc,argv,master,worker_init,worker_done) BIND(C,NAME="mdlLaunch")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR, C_FUNPTR
      INTEGER(KIND=C_INT), VALUE :: argc
      TYPE(C_PTR), VALUE :: argv
      type(c_funptr), intent(in), value :: master
      type(c_funptr), intent(in), value :: worker_init, worker_done
    end subroutine mdl_launch

    subroutine mdl_abort(mdl) BIND(C,NAME="mdlAbort")
      import :: mdl_t
      type(mdl_t)::mdl
    end subroutine mdl_abort

    subroutine mdl_add_service(mdl,sid,p1,service,nin,nout) BIND(C,NAME="mdlAddService")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_FUNPTR
      import :: mdl_t
      type(mdl_t)::mdl
      INTEGER(C_INT), VALUE :: sid
      type(*) :: p1
      type(c_funptr), intent(in), value :: service
!      procedure(mdl_service),pointer,intent(in)::service
      integer(c_int), VALUE :: nin, nout
    end subroutine mdl_add_service

    subroutine mdl_commit_services(mdl) BIND(C,NAME="mdlCommitServices")
      import :: mdl_t
      type(mdl_t)::mdl
    end subroutine mdl_commit_services

    subroutine mdl_handler(mdl) BIND(C,NAME="mdlHandler")
      import :: mdl_t
      type(mdl_t)::mdl
    end subroutine mdl_handler

    integer function mdl_req_service(mdl,who,what,msg,n) BIND(C,NAME="mdlReqService")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      import :: mdl_t
      type(mdl_t)::mdl
      type(C_PTR),VALUE :: msg
      integer(c_int), VALUE :: who, what, n
    end function mdl_req_service

    subroutine mdlGetReply(mdl,rid,msg,n) BIND(C,NAME="mdlGetReply")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT, C_PTR
      import :: mdl_t
      type(mdl_t)::mdl
      integer(c_int), VALUE :: rid
      type(c_ptr), VALUE :: msg
      integer(c_int) :: n
    end subroutine mdlGetReply

    integer function mdl_proc_to_thread(mdl,iProc) BIND(C,NAME="mdlBaseProcToThread")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
      import :: mdl_t
      type(mdl_t)::mdl
      integer(c_int), VALUE :: iProc
    end function mdl_proc_to_thread

    integer function mdl_thread_to_proc(mdl,iThread) BIND(C,NAME="mdlBaseThreadToProc")
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_INT
      import :: mdl_t
      type(mdl_t)::mdl
      integer(c_int), VALUE :: iThread
    end function mdl_thread_to_proc

  end interface

  interface mdl_send_request
    module procedure mdl_send_request_array, mdl_send_request_scalar
  end interface mdl_send_request

  contains
    integer function mdl_threads(mdl)
      type(mdl_t)::mdl
      mdl_threads = mdl%ncpu
    end function mdl_threads

    integer function mdl_self(mdl)
      type(mdl_t)::mdl
      mdl_self = mdl%myid+1 ! CAREFUL: Fortran is 1 based because it is stupid
    end function mdl_self

    integer function mdl_core(mdl)
      type(mdl_t)::mdl
      mdl_core = mdl%icore+1 ! CAREFUL: Fortran is 1 based because it is stupid
    end function mdl_core

    integer function mdl_cores(mdl)
      type(mdl_t)::mdl
      mdl_cores = mdl%ncore
    end function mdl_cores

    integer function mdl_proc(mdl)
      type(mdl_t)::mdl
      mdl_proc = mdl%irank+1 ! CAREFUL: Fortran is 1 based because it is stupid
    end function mdl_proc

    integer function mdl_procs(mdl)
      type(mdl_t)::mdl
      mdl_procs = mdl%nrank
    end function mdl_procs

    integer function mdl_send_request_array(mdl,mdl_function_id,target_cpu,input_size,output_size,input_array)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC
      implicit none
      type(mdl_t)::mdl
      integer,intent(in)::mdl_function_id
      integer,intent(in)::target_cpu,input_size,output_size
      integer,intent(in),dimension(1:input_size),target::input_array

      mdl_send_request_array = mdl_req_service(mdl,target_cpu-1,mdl_function_id,C_LOC(input_array),input_size)
    end function mdl_send_request_array

    integer function mdl_send_request_scalar(mdl,mdl_function_id,target_cpu,input_size,output_size,input)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC
      implicit none
      type(mdl_t)::mdl
      integer,intent(in)::mdl_function_id
      integer,intent(in),optional::target_cpu,input_size,output_size
      TYPE(*),intent(in),optional,target::input
      integer::size
      if (present(input_size)) then
        size = input_size
      else
        size = 0
      end if
      mdl_send_request_scalar = mdl_req_service(mdl,target_cpu-1,mdl_function_id,C_LOC(input),size)
    end function mdl_send_request_scalar

    subroutine mdl_get_reply_scalar(mdl,rid,output_size,output)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_INT
      type(mdl_t)::mdl
      integer(c_int), VALUE :: rid
      type(*), optional, target :: output
      integer(c_int), optional :: output_size
      integer(c_int)::dummy
      call mdlGetReply(mdl,rid,C_LOC(output),dummy)
    end subroutine mdl_get_reply_scalar

    ! This is so totally broken it isn't funny, but we can LINK at least. The version above is "correct-ish"
    subroutine mdl_get_reply_array(mdl,rid,output_size,output_array)
      USE, INTRINSIC :: ISO_C_BINDING, ONLY : C_LOC, C_INT
      implicit none
      type(mdl_t)::mdl
      integer(c_int), VALUE :: rid
      integer::output_size
      integer,dimension(1:output_size),target::output_array
      integer(c_int)::dummy
      call mdlGetReply(mdl,rid,C_LOC(output_array),dummy)
    end subroutine mdl_get_reply_array

end module mdl_module
