pro comp_test, file, Pass=Pass, NoExit=NoExit

Catch, Error_status

if Error_status NE 0 then begin
	print,'Error '+strtrim(string(Error_status),2)+': '+!Error_State.Msg
	Pass = 0
	if NOT Keyword_Set(NoExit) then $
		exit, status=abs(Error_status)
endif	

RESOLVE_ROUTINE, file, /Either, /Compile_Full_File
;RESOLVE_ALL, /Quiet

Pass = 1

end