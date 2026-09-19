      IntegerfunctionLENGTH(Needle)
      Character*(*)Needle
      LENGTH=LEN(Needle)
09999 IF(.NOT.(LENGTH>0.AND.Needle(LENGTH:LENGTH)==' '))GOTO09998
      LENGTH=LENGTH-1
      GOTO09999
09998 CONTINUE
      Return
      end
