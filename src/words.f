      Subroutinewords(string,wordm,lwordm,nwordm)
      character*(*)string,wordm(*)
      integerlwordm(*),nwordm
      character*132dmword
      integeri,lm,lstring
      lstring=len_trim(string)
09999 IF(.NOT.(lstring>0.AND.string(max(lstring,1):max(lstring,1))==' ')
     *)GOTO09998
      lstring=lstring-1
      GOTO09999
09998 CONTINUE
      i=1
      nwordm=0
09996 IF(.NOT.(i<=lstring))GOTO09995
09993 IF(.NOT.(i<=lstring.AND.string(i:i)==' '))GOTO09992
      i=i+1
      GOTO09993
09992 CONTINUE
      lm=0
09990 IF(.NOT.(i<=lstring.AND.string(i:i)/=' '.AND.lm<132))GOTO09989
      lm=lm+1
      dmword(lm:lm)=string(i:i)
      i=i+1
      GOTO09990
09989 CONTINUE
      if(lm>0)then
      nwordm=nwordm+1
      wordm(nwordm)(1:lm)=dmword(1:lm)
      lwordm(nwordm)=lm
      endif
      GOTO09996
09995 CONTINUE
      return
      end
