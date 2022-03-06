! filename: readXsection_3zones_NaturalXsection.f90
! update 20210526
subroutine readXsection(k,lftBnkMann,rmanning_main,rgtBnkMann,leftBnkX_given,rghtBnkX_given,timesDepth,num_reach)

	! k = sequence of x-section from the upstream of a given reach
	! lftBnkMann = Manning's coefficients for the left bank
	! rmanning_main = Manning's coefficients for the main channel
	! rgtBnkMann = Manning's coefficients for the left bank
	! leftBnkX_given = The destance X, at which the channel is divided between left bank and main channel
	! rghtBnkX_given = The destance X, at which the channel is divided between main channel and right bank
	! timesDepth = The multipliyer of the depth upto which level the attribute table will be created
	! num_reach = the sequence of the river reach within a network

    use constants_module
    use arrays_module
    use arrays_section_module
    use xsec_attribute_module

    implicit none
    save

    integer, intent(in) :: k, num_reach
    real, intent(in) :: rmanning_main,lftBnkMann,rgtBnkMann,leftBnkX_given,rghtBnkX_given, timesDepth

    real xcs(maxTableLength), ycs(maxTableLength), el1(nel,3),a1(nel,3),peri1(nel,3),redi1(nel,3),redi1All(nel)
    real conv1(nel,3), tpW1(nel,3), diffArea(nel,3), newI1(nel,3), diffPere(nel,3), newdPdA(nel), diffAreaAll(nel), diffPereAll(nel)
    integer i_start(nel), i_end(nel), i_area, i_find, i, j, jj, num
    real el_min, el_max, el_range, el_incr, el_now, x1, y1, x2, y2, x_start, x_end, waterElev, leftBnkX,rghtBnkX
    real f2m, cal_area, cal_peri, cal_topW, cal_dist, cal_tri_area, cal_multi_area, cal_perimeter, diffAreaCenter
    real compoundSKK(nel), compoundMann
    integer i1, i2

    real allXcs(maxTableLength,3), allYcs(maxTableLength,3)
    real leftBnkY, rghtBnkY,rmanning
    integer mainChanStrt, mainChanEnd, totalNodes(3), kkk, startFound, endFound
    !character*4 file_num

    ! Assign some parameters
    ! f2m is conversion if the data is in feet
    f2m=1.0

!     Open data file
    leftBnkX=leftBnkX_given
    rghtBnkX=rghtBnkX_given

    write(file_num,'(i4.4)')k
      startFound = 0
      endFound = 0



        open(11,file=trim(xSection_path(num_reach))//file_num)
!     Skip headings

      read(11,*)

!     Get CS data
      do i=2,maxTableLength
          read(11,*,end=300)x1,y1
          xcs(i)=x1*f2m
          ycs(i)=y1*f2m
          if ((xcs(i) .ge. leftBnkX) .and. (startFound .eq. 0)) then
            mainChanStrt = i-1
            startFound = 1

          end if
          if ((xcs(i) .ge. rghtBnkX) .and. (endFound .eq. 0)) then
            mainChanEnd = i-1
            endFound = 1
          end if
      enddo
300   close(11)
      num=i
      !print*, leftBnkX, rghtBnkX
      if (leftBnkX .lt. minval(xcs(2:num-1))) leftBnkX = minval(xcs(2:num-1))
      if (rghtBnkX .gt. maxval(xcs(2:num-1))) rghtBnkX = maxval(xcs(2:num-1))
      !print*, leftBnkX, rghtBnkX
      !pause

      leftBnkY = ycs(mainChanStrt)+(leftBnkX-xcs(mainChanStrt))/&
      (xcs(mainChanStrt+1)-xcs(mainChanStrt))*(ycs(mainChanStrt+1)-ycs(mainChanStrt))
      rghtBnkY = ycs(mainChanEnd)+(rghtBnkX-xcs(mainChanEnd))/&
      (xcs(mainChanEnd+1)-xcs(mainChanEnd))*(ycs(mainChanEnd+1)-ycs(mainChanEnd))

      !print*, mainChanStrt, xcs(mainChanStrt+1)
      !print*, mainChanEnd, xcs(mainChanEnd+1
      !print*, leftBnkX, rghtBnkX,leftBnkY, rghtBnkY




!     get max elevation, min elevation, elevation range
      el_min=99999.
      el_max=-99999.
      do i=2,num-1
          if(ycs(i).lt.el_min)el_min=ycs(i)
          if(ycs(i).gt.el_max)el_max=ycs(i)
      enddo
      el_range=(el_max-el_min)*timesDepth
      el_incr=el_range/real(nel-1.0)


      ! re-read cross section data
      open(11,file=trim(xSection_path(num_reach))//file_num)
      !     Skip headings
      read(11,*)
      do i=2,mainChanStrt
          read(11,*,end=400)x1,y1
          allXcs(i,1)=x1*f2m
          allYcs(i,1)=y1*f2m
      enddo

      do i=3,mainChanEnd-mainChanStrt+2
          read(11,*,end=400)x1,y1
          allXcs(i,2)=x1*f2m
          allYcs(i,2)=y1*f2m
      enddo

      do i=3,num
          read(11,*,end=400)x1,y1
          allXcs(i,3)=x1*f2m
          allYcs(i,3)=y1*f2m
      enddo

400   close(11)


      allXcs(mainChanStrt+1,1)=leftBnkX
      allYcs(mainChanStrt+1,1)=leftBnkY

      allXcs(1,1)=allXcs(2,1)
      allYcs(1,1)=el_min+el_range+1.

      allXcs(mainChanStrt+2,1)=allXcs(mainChanStrt+1,1)
      allYcs(mainChanStrt+2,1)=el_min+el_range+1.

      allXcs(2,2)=leftBnkX
      allYcs(2,2)=leftBnkY
      allXcs(mainChanEnd-mainChanStrt+3,2)=rghtBnkX
      allYcs(mainChanEnd-mainChanStrt+3,2)=rghtBnkY

      allXcs(1,2)=allXcs(2,2)
      allYcs(1,2)=el_min+el_range+1.

      allXcs(mainChanEnd-mainChanStrt+4,2)=allXcs(mainChanEnd-mainChanStrt+3,2)
      allYcs(mainChanEnd-mainChanStrt+4,2)=el_min+el_range+1.

      allXcs(2,3) = rghtBnkX
      allYcs(2,3) = rghtBnkY

      allXcs(1,3) = allXcs(2,3)
      allYcs(1,3) = el_min+el_range+1.

      allXcs(i,3) = allXcs(i-1,3)
      allYcs(i,3) = el_min+el_range+1.

      totalNodes(1) = mainChanStrt+2
      totalNodes(2) = mainChanEnd-mainChanStrt+4
      totalNodes(3) = i

      !print*, 'xcs', xcs
      !print*, 'ycs', ycs
      !print*, 'allXcs1', allXcs(1:totalNodes(1),1)
      !print*, 'allYcs1', allYcs(1:totalNodes(1),1)
      !print*, 'allXcs2', allXcs(1:totalNodes(2),2)
      !print*, 'allYcs2', allYcs(1:totalNodes(2),2)
      !print*, 'allXcs3', allXcs(1:totalNodes(3),3)
      !print*, 'allYcs3', allYcs(1:totalNodes(3),3)



      ! new part 20210526
      if (xcs(2) .ge. leftBnkX_given) then
        allXcs(:,1) = 0.
        allYcs(:,1) = 0.
      end if

      if (xcs(num-1) .le. rghtBnkX_given) then
        allXcs(:,3) = 0.
        allYcs(:,3) = 0.
      end if


!     output cs data
      !open(11,file=trim(xSection_path(num_reach))//file_num//'.dat')
      !do i=1,num
      !    write(11,*)xcs(i),ycs(i) !changing all into m unit
      !enddo
      !close(11)

      open(11,file=trim(xSection_path(num_reach))//file_num//'_lines')
      open(22,file=trim(xSection_path(num_reach))//file_num//'_tab')
      write(22,'(140a)')'Elev(m)   Area(m2)   Peri(m)   Radi(m)   Conv(m3/s)   topWidth(m)   newI1(m3)   dPdA(1/m)   rmanning'  ! Hu changed

    xcs = 0.
    ycs = 0.
    newI1=0.0 !Hu changed
    do kkk=1,3
      num = totalNodes(kkk)
      xcs(1:num) = allXcs(1:num,kkk)
      ycs(1:num) = allYcs(1:num,kkk)
      if (kkk .eq. 1) rmanning = lftBnkMann
      if (kkk .eq. 2) rmanning = rmanning_main
      if (kkk .eq. 3) rmanning = rgtBnkMann
      do j=1,nel
          el_now=el_min+real(j-1)*el_incr

          if(abs(el_now - el_min) < TOLERANCE) then
            el_now=el_now+0.00001
          end if

          i_start(1)=-999
          i_end(1)=-999
          i_area=0
          i_find=0
        do i=1,num-1
          y1=ycs(i)
          y2=ycs(i+1)
          if(el_now.le.y1 .and. el_now.gt.y2 .and. i_find.eq.0)then
              i_find=1
              i_area=i_area+1
              i_start(i_area)=i
          endif
          if(el_now.gt.y1 .and. el_now.le.y2 .and. i_find.eq.1)then
              i_find=0
              i_end(i_area)=i
          endif
        enddo

        cal_area=0.
        cal_peri=0.
        cal_topW=0.


        do i=1,i_area
            x1=xcs(i_start(i))
            x2=xcs(i_start(i)+1)
            y1=ycs(i_start(i))
            y2=ycs(i_start(i)+1)
            if(y1.eq.y2)then
                x_start=x1
            else
                x_start=x1+(el_now-y1)/(y2-y1)*(x2-x1)
            endif
            write(11,*)x_start,el_now

            x1=xcs(i_end(i))
            x2=xcs(i_end(i)+1)
            y1=ycs(i_end(i))
            y2=ycs(i_end(i)+1)

            if(y1.eq.y2)then
              x_end=x1
            else
              x_end=x1+(el_now-y1)/(y2-y1)*(x2-x1)
            endif

            cal_topW=x_end-x_start+cal_topW

            write(11,*)x_end,el_now
            write(11,*)'NaN NaN'

            i1=i_start(i)
            i2=i_end(i)

            cal_area = cal_area    &
                     +cal_tri_area(el_now,x_start,xcs(i1+1),ycs(i1+1))    &
                     +cal_multi_area(el_now,xcs,ycs,maxTableLength,i1+1,i2)    &
                     +cal_tri_area(el_now,x_end,xcs(i2),ycs(i2))
            cal_peri = cal_peri    &
                    +cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))    &
                    +cal_perimeter(xcs,ycs,maxTableLength,i1+1,i2)    &
                    +cal_dist(x_end,el_now,xcs(i2),ycs(i2))
            if(i1.eq.1)cal_peri=cal_peri    &
                     -cal_dist(x_start,el_now,xcs(i1+1),ycs(i1+1))
            if(i2.eq.(num-1))cal_peri=cal_peri    &
                     -cal_dist(x_end,el_now,xcs(i2),ycs(i2))

        end do

        el1(j,kkk)=el_now
        a1(j,kkk)=cal_area
        peri1(j,kkk)=cal_peri
        redi1(j,kkk)=a1(j,kkk)/peri1(j,kkk)

        conv1(j,kkk)=1./rmanning*a1(j,kkk)*(redi1(j,kkk))**(2./3.)
        !print*, j, kkk, peri1(j,kkk)
        !print*, j, kkk, conv1(j,kkk)
        if (peri1(j,kkk) .le. TOLERANCE) then
            redi1(j,kkk) =0.0; conv1(j,kkk)=0.0
        endif
        tpW1(j,kkk)=cal_topW

        !if(j.eq.1) then
        if (el_now .le. minval(ycs(1:num))) then
          diffArea(j,kkk)=a1(j,kkk)
          diffPere(j,kkk)=peri1(j,kkk)
        else
          diffArea(j,kkk)=a1(j,kkk)-a1(j-1,kkk)
          diffPere(j,kkk)=peri1(j,kkk)-peri1(j-1,kkk)
        endif

        waterElev=el1(j,kkk)
        do jj=2,j
          diffAreaCenter=el1(jj,kkk)-el_incr*0.5
          !print*,waterElev,diffAreaCenter; pause
          newI1(j,kkk)=newI1(j,kkk)+diffArea(jj,kkk)*(waterElev-diffAreaCenter)
        enddo
        !print*, j, kkk, newI1(j,kkk)
        !if (diffArea(j,kkk) .eq. 0) diffArea(j,kkk) = TOLERANCE
        !newdPdA(j,kkk)=diffPere(j,kkk)/diffArea(j,kkk)



        !print*, j, kkk, newI1(j,kkk)

        !write(22,10)el1(j,kkk),a1(j,kkk),peri1(j,kkk),redi1(j,kkk),conv1(j,kkk),    &
        !           tpW1(j,kkk),newI1(j,kkk),newdPdA(j,kkk)

        end do

      end do

      do j = 1,nel
        el_now=el1(j,1)

        !newdPdA(j)=diffAreaAll(j)/diffPereAll(j)
        if (j .eq. 1) then
            newdPdA(j) = sum(peri1(j,:)) / sum(a1(j,:))
        else
            newdPdA(j)= (sum(peri1(j,:)) - sum(peri1(j-1,:))) / (sum(a1(j,:)) - sum(a1(j-1,:)))
        end if

        !compoundMann = sqrt((abs(tpW1(j,1))*lftBnkMann ** 2. + abs(tpW1(j,2))*rmanning_main ** 2.+&
        ! abs(tpW1(j,3))*rgtBnkMann ** 2.) / (abs(tpW1(j,1))+abs(tpW1(j,2))+abs(tpW1(j,3))))

        compoundMann = sqrt((abs(peri1(j,1))*lftBnkMann ** 2. + abs(peri1(j,2))*rmanning_main ** 2.+&
         abs(peri1(j,3))*rgtBnkMann ** 2.) / (abs(peri1(j,1))+abs(peri1(j,2))+abs(peri1(j,3))))
        compoundSKK(j) = 1. / compoundMann

        !print*, abs(tpW1(j,1)), abs(tpW1(j,2)), abs(tpW1(j,3)),lftBnkMann,rmanning_main,rgtBnkMann
        !pause

        redi1All(j)=sum(a1(j,:)) /sum(peri1(j,:))






        !print*, sum(newI1(j,:))
        !write(22,10)el1(j,1),sum(a1(j,:)),sum(peri1(j,:)),redi1All(j),sum(conv1(j,:)),    &
        !   sum(tpW1(j,:)),sum(newI1(j,:)),newdPdA(j)
        write(22,10)el1(j,1),sum(a1(j,:)),sum(peri1(j,:)),redi1All(j),sum(conv1(j,:)),    &
           abs(tpW1(j,1))+abs(tpW1(j,2))+abs(tpW1(j,3)),sum(newI1(j,:)),newdPdA(j), compoundSKK(j)

10      format(f9.2,3f12.2,2f20.3,3f16.3)
        ! new change 20200107 to make the model faster
        xsec_tab(1,j,k,num_reach) = el1(j,1)
        xsec_tab(2,j,k,num_reach) = sum(a1(j,:))
        xsec_tab(3,j,k,num_reach) = sum(peri1(j,:))
        xsec_tab(4,j,k,num_reach) = redi1All(j)
        xsec_tab(5,j,k,num_reach) = sum(conv1(j,:))
        xsec_tab(6,j,k,num_reach) = abs(tpW1(j,1))+abs(tpW1(j,2))+abs(tpW1(j,3))
        xsec_tab(7,j,k,num_reach) = sum(newI1(j,:))
        xsec_tab(8,j,k,num_reach) = newdPdA(j)
        xsec_tab(11,j,k,num_reach) = compoundSKK(j)
        !print*,  el1(j,1), a1(j,1)+a1(j,2)+a1(j,3), peri1(j,1)+peri1(j,2)+peri1(j,3)
        !pause

     ! print*,j,peri1(j,1),peri1(j,2),peri1(j,3), a1(j,1),a1(j,2),a1(j,3)


      end do
!pause 2030

      close(11)
      close(22)


      z(k,num_reach)= el_min
      !pause 3030
end subroutine


function cal_tri_area(el,x0,x1,y1)
      cal_tri_area=abs(0.5*(x1-x0)*(el-y1))
      return
end function

function cal_trap_area(el,x1,y1,x2,y2)
      cal_trap_area=abs(0.5*(x2-x1)*(el-y1+el-y2))
      return
end function

function cal_multi_area(el,xx,yy,n,i1,i2)
      real xx(n),yy(n)

      area=0
      do i=i1,i2-1
          x1=xx(i)
          y1=yy(i)
          x2=xx(i+1)
          y2=yy(i+1)
          area=area+cal_trap_area(el,x1,y1,x2,y2)
      enddo
      cal_multi_area=area
      return
end function

function cal_dist(x1,y1,x2,y2)
          cal_dist=sqrt((x1-x2)*(x1-x2)+(y1-y2)*(y1-y2)+1.e-32)
      return
end function

function cal_perimeter(xx,yy,n,i1,i2)
        real xx(n),yy(n)

      p=0.
      do i=i1,i2-1
        x1=xx(i)
        y1=yy(i)
        x2=xx(i+1)
        y2=yy(i+1)
        p=p+cal_dist(x1,y1,x2,y2)
      enddo
      cal_perimeter=p
      return
end function
