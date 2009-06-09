pro link_test, file, target, max_agrs, min_args, Pass=Pass, NoExit=NoExit

Catch, Error_status

if Error_status NE 0 then begin
	print,'Error '+strtrim(string(Error_status),2)+': '+!Error_State.Msg
	Pass = 0
	if NOT Keyword_Set(NoExit) then $
		exit, status=255
endif

linkimage, target, file, 1, target, max_args=max_args, min_args=min_args

Pass = 1

end
