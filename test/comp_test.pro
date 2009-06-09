pro comp_test, file, Pass=Pass, NoExit=NoExit

Catch, Error_status

if Error_status NE 0 then begin
	print,'Error '+strtrim(string(Error_status),2)+': '+!Error_State.Msg
	Pass = 0
	if NOT Keyword_Set(NoExit) then $
		exit, status=255
endif	

RESOLVE_ROUTINE, file
;RESOLVE_ALL, /Quiet

Pass = 1

end