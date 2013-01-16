!! Module for unified error messaging
!! -
!! Fortran 90/95 lacks proper error handling.
!! This module tries to overcome with that problem.
!! Clients (program modules) that use this module can
!! register itselves with unique error id, and then
!! issue errors using that id. Upper-level program blocks
!! can check the error status and choose wether to stop
!! or continue, and where to print the error message.
!! -
!! Error messages include the id of the issuer, integer
!! error code and a string message.
!! -
!! If you use this code in a publication, please make a reference to:
!! A. Penttilä, Fortran 95 error handling implementation (computer code),
!! http://wiki.helsinki.fi/display/~aipentti@helsinki.fi/Collection+of+codes (2012).
!! -
!! Antti Penttilä
!! 2012
!! Department of Physics, University of Helsinki
MODULE ErrorCapsule

	IMPLICIT NONE
	
	PRIVATE

  !! Private, maximum number of clients
	INTEGER, PARAMETER :: max_clients = 100
  !! Public, maximum message length
  INTEGER, PARAMETER :: max_msg_length = 512
  INTEGER, PARAMETER :: error_too_much_clients = -121
	
	INTEGER :: i, j, k, client_i = 2
	INTEGER, DIMENSION(max_clients) :: error_codes = 0
	
	CHARACTER(LEN=1024) :: str
	CHARACTER(LEN=max_msg_length), DIMENSION(max_clients) :: error_msgs, &
		client_names = "ErrorCapsule"

	INTEGER :: last_error = 0

	LOGICAL :: error_state = .FALSE.

	PUBLIC :: error_state, last_error, get_err_id, issue_error, clear_errors, get_error, &
		print_last_error, print_errors, check_errors, max_msg_length
	
CONTAINS

! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Client must get a client_id first.
FUNCTION get_err_id(name) RESULT(id)
	IMPLICIT NONE
  !! Input, name of the calling program unit/module
	CHARACTER(LEN=*), INTENT(IN) :: name
  !! Result, unique id number for the client
	INTEGER :: id

	IF(client_i >= max_clients) THEN
		id = 0
		CALL issue_error(1,error_too_much_clients, "Too much clients")
	ELSE
		client_i = client_i+1
		id = client_i
		client_names(client_i) = name
	END IF
	
END FUNCTION get_err_id


! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Issue an error using client id
SUBROUTINE issue_error(id, code, msg)
	IMPLICIT NONE
  !! Inputs, client id and integer error code.
	INTEGER, INTENT(IN) :: id, code
  !! Input, string error message for the user.
	CHARACTER(LEN=*), INTENT(IN) :: msg

	last_error = id
	error_codes(id) = code
	error_msgs(id) = msg
	error_state = .TRUE.
  
END SUBROUTINE issue_error


! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Clear errors, use when the possible errors
!! have been handled somehow.
SUBROUTINE clear_errors()
	IMPLICIT NONE

	last_error = 0
	error_codes = 0
	error_msgs = ""
	error_state = .FALSE.
  
END SUBROUTINE clear_errors


! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Get (last) error from specific client id
SUBROUTINE get_error(id, code, msg)
	IMPLICIT NONE
  !! Input, client id
	INTEGER, INTENT(IN) :: id
  !! Output, error code
	INTEGER, INTENT(INOUT) :: code
  !! Output, error message
	CHARACTER(LEN=max_msg_length), INTENT(INOUT) :: msg

	code = error_codes(id)
	msg = error_msgs(id)

END SUBROUTINE get_error


! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Prints last error. Can print to
!! std_out or an open file channel
SUBROUTINE print_last_error(ch_id)
	IMPLICIT NONE
  !! Input, use 0 to print to std_out, and 
  !! positive integer to print to that unit.
  !! Be sure that the unit is open for writing.
	INTEGER, INTENT(IN) :: ch_id

	IF(ch_id > 0) THEN
		WRITE(ch_id,*) "***Last error from ErrorCapsule:***"
	ELSE
		WRITE(*,*) "***Last error from ErrorCapsule:***"
	END IF
	WRITE(str, *) "[FROM ", TRIM(client_names(last_error)), &
		"]: [code ", error_codes(last_error), "]: ", TRIM(error_msgs(last_error))
	IF(ch_id > 0) THEN
		WRITE(ch_id,*) TRIM(str)
	ELSE
		WRITE(*,*) TRIM(str)
	END IF
	
END SUBROUTINE print_last_error


! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Prints all active errors, i.e. errors
!! that have not been cleared. Prints to
!! std_out or an open file channel
SUBROUTINE print_errors(ch_id)
	IMPLICIT NONE
  !! Input, use 0 to print to std_out, and 
  !! positive integer to print to that unit.
  !! Be sure that the unit is open for writing.
	INTEGER, INTENT(IN) :: ch_id

	IF(ch_id > 0) THEN
		WRITE(ch_id,*) "***Error states from ErrorCapsule:***"
	ELSE
		WRITE(*,*) "***Error states from ErrorCapsule:***"
	END IF
	DO i=1,client_i
		IF(error_codes(i) /= 0) THEN
			WRITE(str, *) "[FROM ", TRIM(client_names(i)), &
				"]: [code ", error_codes(i), "]: ", TRIM(error_msgs(i))
			IF(ch_id > 0) THEN
				WRITE(ch_id,*) TRIM(str)
			ELSE
				WRITE(*,*) TRIM(str)
			END IF
		END IF
	END DO
	
END SUBROUTINE print_errors


! PUBLIC!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Check and possibly print out error(s).
!! Optionally can stop the execution of the program if
!! there have been errors (default). Can print errors to
!! std_out (default) or to an open file channel.
SUBROUTINE check_errors(END_EXEC, CHANNEL)
	IMPLICIT NONE
  !! Optional input. If .TRUE. (default), will stop the execution
  !! of the whole program if errors are found.
	LOGICAL, OPTIONAL :: END_EXEC
  !! Optional nput, use 0 to print to std_out, and 
  !! positive integer to print to that unit.
  !! Be sure that the unit is open for writing.
  INTEGER, INTENT(IN), OPTIONAL :: CHANNEL
  INTEGER :: ch_id
	LOGICAL :: full_stop
	
	IF(PRESENT(END_EXEC)) THEN
		full_stop = END_EXEC
	ELSE
		full_stop = .TRUE.
	END IF
  
  IF(PRESENT(CHANNEL)) THEN
    ch_id = CHANNEL
  ELSE
    ch_id = 0
  END IF

	IF(error_state) THEN
		error_state = .FALSE.
		CALL print_errors(ch_id)
		IF(full_stop) STOP
	END IF
	
END SUBROUTINE check_errors


END MODULE ErrorCapsule