! Utilities to read information from files and store them for use by the
! main program. 
! This file contains the following modules:
! 1. IniFile   -- Reads in name/value pairs from a parameter file, with 
!                 each line of the form 'name = array of values'.
!                 Adapted from the camb parser at http://camb.info
! 2. Data_Repo -- Reads in tabular data from files and stores it in a 
!                 repository. It assumes each row is of the form
!                 'indep_var dep_var(1) dep_var(2) ... dep_var(n)'
! and the following subroutines:
! 1. writearray -- outputs arrays of independent and dependent 
!                  variables to a file/stdout
!---------------------------------------------------------------------
      subroutine writearray(x, y, m, n, xstride, ystride, header, fp)!{{{
      !*************************************************************
      ! Subroutine to write to a file indexed by fp the values 
      ! x(xstride+1:xstride+m) (abscissae), and y(ystride+1:ystride+m*n) 
      ! (n columns of ordinates) in double precision. Works for n-d arrays 
      ! stored in column major. This form makes it easy to output everything. 
      ! Typecast inputs to this function, if so needed!
        
      Integer m, n, xstride, ystride, fp!{{{
      Double precision x(xstride+m), y(ystride + m*n)
      Character*(*) header
      
      Integer i, j!}}}

      Write(fp,*) '#',header

      i=1
      j=1

      ! Gymnastics to avoid spurious newline at end.
      ! First write upto penultimate row
      Do while (i.le.m-1)
         Write(fp,100) x(xstride+i)
         ! Write upto penultimate element, array stored in column major.
         Do while (j.le.n-1)
            Write(fp,100) y(ystride + (j-1)*m + i)
            j = j + 1
         End do
         ! Write the last element
         Write(fp,110) y(ystride + (n-1)*m + i)
         j = 1
         i = i + 1
      End do
      ! Now write the last row
      Write(fp,100) x(xstride+m)
      Do j=1, n
         Write(fp,100) y(ystride + j*m)
      End do

 100  Format(1pE24.16, $)
 110  Format(1pE24.16)

      Return
      !*************************************************************
      End Subroutine ! writearray!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
module IniFile!{{{
  Use defaults
  implicit none
  public

  ! Blocks to increase sizes in reading ini file
  integer, parameter :: Ini_blocksize = 128

  type TNameValue!{{{
    ! Each entry is possibly an array
    integer nvals
    integer Capacity
    character(INI_MAX_NAME_LEN)  :: Name
    character(INI_MAX_STRING_LEN), dimension(:), pointer :: Values
  end type TNameValue!}}}

  Type TNameValueList!{{{
    integer Count
    integer Capacity
    type(TNameValue), dimension(:), pointer :: Items
  end Type TNameValueList!}}}

  Type TIniFile!{{{
   logical HashComments
   Type (TNameValueList) :: L, ReadValues
  end Type TIniFile!}}}
 
  Type(TIniFile) :: DefIni

contains

   subroutine TNameValue_Init(L)!{{{
    Type (TNameValue) :: L
    
     L%nvals = 0
     L%Capacity = 0
     L%Name = ''
     nullify(L%Values)

   end subroutine TNameValue_Init!}}}

   subroutine TNameValueList_Init(L)!{{{
    Type (TNameValueList) :: L
    
     L%Count = 0
     L%Capacity = 0
     nullify(L%Items)

   end subroutine TNameValueList_Init!}}}

   subroutine TNameValue_Clear(L)!{{{
    Type (TNameValue) :: L
    integer status

    deallocate (L%Values, stat = status)
    call TNameValue_Init(L)

   end subroutine TNameValue_Clear!}}}

   subroutine TNameValueList_Clear(L)!{{{
    Type (TNameValueList) :: L
    integer i, status
     
    do i=L%count,1,-1
     call TNameValue_Clear(L%Items(i))
    end do
    deallocate (L%Items, stat = status)
    call TNameValueList_Init(L)

   end subroutine TNameValueList_Clear!}}}

   subroutine TNameValueList_ValueOf(L, AName, index, AValue)!{{{
   ! Possibly asking for value in array
   ! If index > size, return blank
     Type (TNameValueList), intent(in) :: L
     character(LEN=*), intent(in) :: AName
     Integer, intent(in) :: index
     Character(LEN=*), intent(out) :: AValue

     integer i

     do i=1, L%Count
       if (L%Items(i)%Name == AName) then
         if ((index<=L%Items(i)%nvals) .and. (index>0)) then
           AValue = L%Items(i)%Values(index)
           return
         end if
       end if
     end do
     AValue = ''

   end subroutine TNameValueList_ValueOf!}}}
   
   function TNameValueList_IndexOf(L, AName) result (AValue)!{{{
     Type (TNameValueList), intent(in) :: L
     character(LEN=*), intent(in) :: AName
     integer :: AValue
     integer i

     do i=1, L%Count
       if (L%Items(i)%Name == AName) then
          AValue = i
          return
       end if
     end do
     AValue = -1
     
   end function TNameValueList_IndexOf!}}}

   subroutine TNameValueList_NumberOf(L, AName, AValue)!{{{
   ! Find number of entries 
     Type (TNameValueList), intent(in) :: L
     character(LEN=*), intent(in) :: AName
     Integer, intent(out) :: AValue

     integer i

     i = TNameValueList_IndexOf(L, Aname)

     If (i>0) then
       AValue = L%Items(i)%nvals
     Else
       AValue = -1
     End If
     return

   end subroutine TNameValueList_NumberOf!}}}

   function TNameValueList_HasKey(L, AName) result (AValue)!{{{
     Type (TNameValueList), intent(in) :: L
     character(LEN=*), intent(in) :: AName
     logical :: AValue
    
     AValue = TNameValueList_IndexOf(L,AName) /= -1
     
   end function TNameValueList_HasKey!}}}
    
   subroutine TNameValueList_Add(L, AName, AValue, only_if_undefined)!{{{
    Type (TNameValueList) :: L
    character(LEN=*), intent(in) :: AName, AValue
    logical, optional, intent(in) :: only_if_undefined

    logical isDefault
    integer :: nameindex
    integer :: valueindex
   
    if (present(only_if_undefined)) then
     isDefault = only_if_undefined
    else
     isDefault = .true.
    end if

    if (TNameValueList_HasKey(L,AName)) then
      ! Merely add to list of entries
      ! Find index for the name
      nameindex = TNameValueList_IndexOf(L, AName)
      if (L%Items(nameindex)%nvals == L%Items(nameindex)%Capacity) then
        call TNameValue_SetCapacity(L%Items(nameindex), L%Items(nameindex)%Capacity + Ini_blocksize)
      end if
      valueindex = L%Items(nameindex)%nvals + 1 ! Adding one more to the values
      L%Items(nameindex)%nvals = valueindex
      L%Items(nameindex)%Values(valueindex) = Avalue
    else ! create new key and add entry
      if (L%Count == L%Capacity) call TNameValueList_SetCapacity(L, L%Capacity + Ini_blocksize)
      L%Count = L%Count + 1
      call TNameValue_Init(L%Items(L%Count)) ! Initialize field
      call TNameValue_SetCapacity(L%Items(L%Count), 1) ! set the capacity to one
      L%Items(L%Count)%nvals =  1 ! starts out with one
      L%Items(L%Count)%Name = AName
      L%Items(L%Count)%Values(1) = AValue
    end if

   end subroutine TNameValueList_Add!}}}

   subroutine TNameValue_SetCapacity(L, C)!{{{
   ! Set size of each entry's array
    Type (TnameValue) :: L
    Integer C
    character(INI_MAX_STRING_LEN), dimension(:), pointer :: TmpValues

    if (L%nvals > 0) then
      if (C < L%nvals) stop 'TnameValue_SetCapacity: smaller than nvals'
      allocate(TmpValues(L%nvals))
      TmpValues = L%Values(1:L%nvals)
      deallocate(L%Values)
      allocate(L%Values(C))
      L%Values(1:L%nvals) = TmpValues
      deallocate(TmpValues)
    else
      allocate(L%Values(C))
    end if
    L%Capacity = C
    
   end subroutine TnameValue_SetCapacity!}}}

   subroutine TNameValueList_SetCapacity(L, C)!{{{
    Type (TNameValueList) :: L
    integer C
    type(TNameValue), dimension(:), pointer :: TmpItems
    
    if (L%Count > 0) then
      if (C < L%Count) stop 'TNameValueList_SetCapacity: smaller than Count'
      allocate(TmpItems(L%Count))
      TmpItems = L%Items(1:L%Count)
      deallocate(L%Items)
      allocate(L%Items(C))
      L%Items(1:L%Count) = TmpItems
      deallocate(TmpItems)
    else
      allocate(L%Items(C))
    end if  
    L%Capacity = C
  
   end subroutine TNameValueList_SetCapacity!}}}

   subroutine TNameValueList_Delete(L, i)!{{{
    Type (TNameValueList) :: L
    integer, intent(in) :: i

    call TNameValue_Clear(L%Items(i))
    if (L%Count > 1) L%Items(i:L%Count-1) = L%Items(i+1:L%Count)
    L%Count = L%Count -1
     
   end subroutine TNameValueList_Delete!}}}

  subroutine Ini_NameValue_Add(Ini,AInLine,only_if_undefined)!{{{
    Type(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: AInLine
    integer EqPos, Hashpos, lastpos, cpos
    logical, optional, intent(in) :: only_if_undefined
    logical isDefault
    character (LEN=len(AInLine)) :: AName, S, InLine

    if (present(only_if_undefined)) then
      isDefault = only_if_undefined  
    else
      isDefault = .false.  
    end if
    InLine=trim(adjustl(AInLine))
    EqPos = scan(InLine,'=')
    if (EqPos==0 .or. InLine(1:1)=='#' .or. InLine(1:7)=='COMMENT' ) return
    AName = trim(adjustl(InLine(1:EqPos-1)))
    if (AName=='') return ! Don't take in empty names
    S = adjustl(InLine(EqPos+1:))

    if (Ini%HashComments) then
      Hashpos=scan(S,'#')
      if (Hashpos /= 0) then
        S  = S(1:Hashpos-1)
      end if
    end if

    lastpos=len_trim(S)
    if (lastpos>1) then
      if (S(1:1)=='''' .and. S(lastpos:lastpos)=='''') then
        S = S(2:lastpos-1)
      end if
    end if
    ! Now scan through line and pass comma separated values
    Do
      cpos = Scan(S,',')
      If (cpos>1) then
        call TNameValueList_Add(Ini%L, AName, trim(adjustl(S(1:cpos-1))),only_if_undefined = isDefault)
        S = S(cpos+1:)
      Else
        call TNameValueList_Add(Ini%L, AName, trim(adjustl(S)),only_if_undefined = isDefault)
        exit
      End if
    End do
    return
  end subroutine Ini_NameValue_Add!}}}

  subroutine Ini_Open(filename, unit_id,  error, Hash_comments)!{{{
     character (LEN=*), intent(IN) :: filename
     integer, intent(IN) :: unit_id
     logical, optional, intent(OUT) :: error
     logical, optional, intent(IN) :: Hash_comments
     logical aerror
          
     if (present(Hash_comments)) then
      call Ini_Open_File(DefIni,filename,unit_id,aerror,Hash_comments)
     else
      call Ini_Open_File(DefIni,filename,unit_id,aerror)
     end if

     if (present(error)) then
       error = aerror
     else
      if (aerror) then
        write (*,*) 'Ini_Open: Error opening file ' // trim(filename)
        stop
      end if
     end if

  end subroutine Ini_Open!}}}

  function Ini_ExtractFilePath(aname)!{{{
    character(LEN=*), intent(IN) :: aname
    character(LEN=INI_MAX_STRING_LEN) Ini_ExtractFilePath
    integer len, i

    len = len_trim(aname)
    do i = len, 1, -1
       if (aname(i:i)=='/') then
          Ini_ExtractFilePath = aname(1:i)
          return
       end if
    end do
    Ini_ExtractFilePath = ''

  end function Ini_ExtractFilePath!}}}

  recursive subroutine Ini_Open_File(Ini, filename, unit_id,            &!{{{
     &                               error, Hash_comments, append,only_if_undefined)
     Type(TIniFile) :: Ini

     character (LEN=*), intent(IN) :: filename
     integer, intent(IN) :: unit_id
     logical, intent(OUT) :: error
     logical, optional, intent(IN) :: Hash_comments
     logical, optional, intent(in) :: append, only_if_undefined

     character (LEN=INI_MAX_STRING_LEN) :: InLine
     integer lastpos, i
     logical doappend, isDefault
     
     if (present(append)) then
      doappend=append
     else
      doappend=.false.
     end if  
     
     if (present(only_if_undefined)) then
        isDefault = only_if_undefined
     else
        isDefault = .false.
     end if
     
     if (.not. doappend) then
       call TNameValueList_Init(Ini%L)
       call TNameValueList_Init(Ini%ReadValues)
     end if
     
    if (present(Hash_comments)) then
     Ini%HashComments = Hash_comments
    else
     Ini%HashComments = .false.
    end if
         
    open(unit=unit_id,file=filename,form='formatted',status='old', err=500)
   
    do 
      read (unit_id,'(a)',end=400) InLine
      if (InLine == 'END') exit;
      if (InLine /= '') then
       call Ini_NameValue_Add(Ini,InLine, only_if_undefined=isDefault) 
      end if
    end do

400 close(unit_id)
    error=.false.
    If (IO_VERBOSE) Print *,"Read in ",Ini%L%Count," fields from ",filename 
    return

500 error=.true.
    return
  end subroutine Ini_Open_File!}}}

  subroutine Ini_Open_Fromlines(Ini, Lines, NumLines, Hash_comments)!{{{
    Type(TIniFile) :: Ini

    integer, intent(IN) :: NumLines
    character (LEN=*), dimension(NumLines), intent(IN) :: Lines
    logical, intent(IN) :: Hash_comments
    integer i

    call TNameValueList_Init(Ini%L)
    call TNameValueList_Init(Ini%ReadValues)

    Ini%HashComments = Hash_comments

    do i=1,NumLines
       call Ini_NameValue_Add(Ini,Lines(i))
    end do

  end  subroutine Ini_Open_Fromlines!}}}

   function Ini_NumberOf(AName) result (AValue)!{{{
   ! Find number of entries 
     character(LEN=*), intent(in) :: AName
     Integer :: AValue

     AValue = Ini_NumberOf_File(DefIni, AName)
   end function Ini_NumberOf!}}}

   function Ini_NumberOf_File(Ini, AName) result (AValue)!{{{
   ! Find number of entries
     Type(TIniFile) :: Ini
     character(LEN=*), intent(in) :: AName
     Integer :: AValue

     call TNameValueList_NumberOf(Ini%L, AName, AValue)
   end function Ini_NumberOf_File!}}}

  subroutine Ini_Close!{{{

    call Ini_close_File(DefIni)

  end subroutine Ini_Close!}}}

  subroutine Ini_Close_File(Ini)!{{{
    Type(TIniFile) :: Ini
   
    call TNameValueList_Clear(Ini%L)
    call TNameValueList_Clear(Ini%ReadValues)

  end  subroutine Ini_Close_File!}}}
  
  function Ini_HasKey(Key) result(AValue)!{{{
   character (LEN=*), intent(IN) :: Key
   logical AValue

   AValue = Ini_HasKey_File(DefIni, Key)

  end function Ini_HasKey!}}}

  function Ini_HasKey_File(Ini, Key) result(AValue)!{{{
   type(TIniFile), intent(in) :: Ini
   character (LEN=*), intent(IN) :: Key
   logical AValue

      Avalue = TNameValueList_HasKey(Ini%L, Key)
      
  end function Ini_HasKey_File!}}}

  function Ini_Read_String(Key, index) result (AValue)!{{{
   character (LEN=*), intent(IN) :: Key
   integer, optional, intent(IN) :: index
   character(LEN=INI_MAX_STRING_LEN) :: AValue

     if (present(index)) then
      AValue = Ini_Read_String_File(DefIni, Key, index)
     else
      AValue = Ini_Read_String_File(DefIni, Key)
     end if

  end function Ini_Read_String!}}}

  function Ini_Read_String_File(Ini, Key, index) result (AValue)!{{{
   Type(TIniFile) :: Ini
   character (LEN=*), intent(IN) :: Key
   Integer, optional, intent(IN) :: index
   character(LEN=INI_MAX_STRING_LEN) :: AValue

   integer :: readindex

   if (present(index)) then
     readindex = index
   else
     readindex = 1 ! default to scalar
   end if

   call TNameValueList_ValueOf(Ini%L, Key, readindex, AValue)

   if (AValue/='') then
     call  TNameValueList_Add(Ini%ReadValues, Key, AValue)
     return
   else if (INI_FAIL_ON_NOT_FOUND) then
     write(*,*) 'key not found : '//trim(Key)
     stop
   else
     return
   end if

  end function Ini_Read_String_File!}}}
  
  function Ini_Read_Int(Key, index) result (AValue)!{{{
   character (LEN=*), intent(IN) :: Key
   integer, optional, intent(IN) :: index
   Integer AValue

     if (present(index)) then
      AValue = Ini_Read_Int_File(DefIni, Key, index)
     else
      AValue = Ini_Read_Int_File(DefIni, Key)
     end if

  end function Ini_Read_Int!}}}

  function Ini_Read_Int_File(Ini, Key, index) result (AValue)!{{{
    Type(TIniFile) :: Ini
    character  (LEN=*), intent(IN) :: Key
    integer, optional, intent(IN) :: index
    integer AValue
 
    character(LEN=INI_MAX_STRING_LEN) :: S
    
    if (present(index)) then
      S = Ini_Read_String_File(Ini, Key, index)
    else
      S = Ini_Read_String_File(Ini, Key)
    end if
 
    call TNameValueList_Add(Ini%ReadValues, Key, S)
    if (verify(trim(S),'-+0123456789') /= 0) goto 10
    read (S,*, err = 10) AValue
    return
10  write (*,*) 'error reading integer for key: '//Key
    stop
  
  end function Ini_Read_Int_File!}}}

  function Ini_Read_Double(Key, index) result (AValue)!{{{
    character (LEN=*), intent(IN) :: Key
    integer, optional, intent(IN) :: index

    double precision AValue

    if (present(index)) then
      AValue = Ini_Read_Double_File(DefIni, Key, index)
    else
      AValue = Ini_Read_Double_File(DefIni, Key)
    end if
  
  end function Ini_Read_Double!}}}

  function Ini_Read_Double_File(Ini,Key,index) result (AValue)!{{{
    Type(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key
    integer, optional, intent(IN) :: index
    Double precision AValue

    character(LEN=INI_MAX_STRING_LEN) :: S

    if (present(index)) then
      S = Ini_Read_String_File(Ini,Key,index)
    else
      S = Ini_Read_String_File(Ini,Key)
    end if

    call TNameValueList_Add(Ini%ReadValues, Key, S)
    read (S,*, err=10) AValue
    return
10  write (*,*) 'error reading double for key: '//Key
    stop

  end function Ini_Read_Double_File!}}}

  function Ini_Read_Real(Key, index) result (AValue)!{{{
    character (LEN=*), intent(IN) :: Key
    integer, optional, intent(IN) :: index

    real AValue

    if (present(index)) then
      AValue = Ini_Read_Real_File(DefIni, Key, index)
    else
      AValue = Ini_Read_Real_File(DefIni, Key)
    end if
  
  end function Ini_Read_Real!}}}

  function Ini_Read_Real_File(Ini,Key,index) result (AValue)!{{{
    Type(TIniFile) :: Ini
    character (LEN=*), intent(IN) :: Key
    integer, optional, intent(IN) :: index
    real AValue

    character(LEN=INI_MAX_STRING_LEN) :: S

    if (present(index)) then
      S = Ini_Read_String_File(Ini,Key,index)
    else
      S = Ini_Read_String_File(Ini,Key)
    end if

    call TNameValueList_Add(Ini%ReadValues, Key, S)
    read (S,*, err=10) AValue
    return
10  write (*,*) 'error reading real for key: '//Key
    stop

  end function Ini_Read_Real_File!}}}

  function Ini_Read_Logical(Key, index) result (AValue)!{{{
    character (LEN=*), intent(IN) :: Key
    integer, optional, intent(IN) :: index

    logical AValue

    if (present(index)) then
      AValue = Ini_Read_Logical_File(DefIni, Key, index)
    else
      AValue = Ini_Read_Logical_File(DefIni, Key)
    end if
  
  end function Ini_Read_Logical!}}}

  function Ini_Read_Logical_File(Ini, Key,index) result (AValue)!{{{
    Type(TIniFile) :: Ini
    character  (LEN=*), intent(IN) :: Key
    integer, optional, intent(IN) :: index
    logical AValue

    character(LEN=INI_MAX_STRING_LEN) :: S

    if (present(index)) then
      S = Ini_Read_String_File(Ini,Key,index)
    else
      S = Ini_Read_String_File(Ini,Key)
    end if

    call  TNameValueList_Add(Ini%ReadValues, Key, S)

    if (verify(trim(S),'10TF') /= 0) goto 10  
    read (S,*, err = 10) AValue
    return

10 write (*,*) 'error reading logical for key: '//Key
   stop
  end function Ini_Read_Logical_File!}}}

  subroutine Ini_SaveReadValues(afile,unit_id)!{{{
   character(LEN=*)  :: afile
   integer, intent(in) :: unit_id

   call Ini_SaveReadValues_File(DefIni, afile, unit_id)

  end subroutine Ini_SaveReadValues!}}}

  subroutine Ini_SaveReadValues_File(Ini, afile, unit_id)!{{{
   Type(TIniFile) :: Ini
   character(LEN=*), intent(in) :: afile
   integer, intent(in) :: unit_id
   integer i, j

   open(unit=unit_id,file=afile,form='formatted',status='replace', err=500)

   do i=1, Ini%ReadValues%Count
     do j=1, Ini%ReadValues%Items(i)%nvals
       write (unit_id,'(a)') trim(Ini%ReadValues%Items(i)%Name) // ' = ' &
                        //trim(Ini%ReadValues%Items(i)%Values(j))
     end do
   end do
   close(unit_id)
   return
500 write(*,*) 'Ini_SaveReadValues_File: Error creating '//trim(afile)

  end subroutine Ini_SaveReadValues_File!}}}

end module IniFile!}}}
!---------------------------------------------------------------------

!---------------------------------------------------------------------
module Data_Repo!{{{
  Use defaults
  implicit none
  public

  ! Blocks to increase data size while reading data file
  integer, parameter :: Tab_blocksize = 128

  type Gen_Array!{{{
    ! A general array with its length      
    integer :: length
    integer :: Capacity
    double precision, dimension(:), pointer :: Values
  end type Gen_Array!}}}

  type Data_Set!{{{
    ! A dataset with an independent variable, a number of dependent 
    ! variables whose values are in a 2d array (such that this matrix 
    ! `printed out' gives the table in the file), and an auxiliary 
    ! array, if needed
    character(INI_MAX_NAME_LEN) :: Name
    integer :: n_indep
    integer :: Capacity
    integer :: n_dep
    double precision, dimension(:), pointer :: Indep_Values
    double precision, dimension(:,:), pointer :: Dep_Values
    type(Gen_Array) :: Aux_Values
  end type Data_Set!}}}

  Type Data_Book!{{{
    ! A collection of datasets
    integer :: Count
    integer :: Capacity
    logical HashComments
    type(Data_Set), dimension(:), pointer :: Datasets
  end Type Data_Book!}}}

  Type(Data_Book) :: Def_Data

  Save
  Contains
    subroutine Data_Repo_Init!{{{
      call Data_Book_Init(Def_Data)
      if (IO_VERBOSE) Print *,"Initialized Data repository"
    end subroutine Data_Repo_Init!}}}

    subroutine Data_Repo_Close!{{{
      call Data_Book_Clear(Def_Data)
    end subroutine Data_Repo_Close!}}}

   subroutine Gen_Array_Init(L)!{{{
    Type (Gen_Array) :: L
     L%length = 0
     L%Capacity = 0
     nullify(L%Values)
   end subroutine Gen_Array_Init!}}}

   subroutine Data_Set_Init(L)!{{{
    Type (Data_Set) :: L
     L%Name = ''
     L%n_indep = 0
     L%Capacity = 0
     L%n_dep = 0
     nullify(L%Indep_Values)
     nullify(L%Dep_Values)
     call Gen_Array_Init(L%Aux_Values)
   end subroutine Data_Set_Init!}}}

   subroutine Data_Book_Init(L)!{{{
    Type (Data_Book) :: L
     L%Count = 0
     L%Capacity = 0
     L%HashComments = .true.
     nullify(L%Datasets)
   end subroutine Data_Book_Init!}}}

   subroutine Gen_Array_Clear(L)!{{{
    Type (Gen_Array) :: L
    integer status
    deallocate (L%Values, stat = status)
    call Gen_Array_Init(L)
   end subroutine Gen_Array_Clear!}}}

   subroutine Data_Set_Clear(L)!{{{
    Type (Data_Set) :: L
    integer status
    deallocate (L%Indep_Values, stat = status)
    deallocate (L%Dep_Values, stat = status)
    call Gen_Array_Clear(L%Aux_Values)
    call Data_Set_Init(L)
   end subroutine Data_Set_Clear!}}}

   subroutine Data_Book_Clear(L)!{{{
    Type (Data_Book) :: L
    integer i, status
    Do i=1,L%Count
      call Data_Set_Clear(L%Datasets(i))
    End do
    deallocate (L%Datasets, stat = status)
    call Data_Book_Init(L)
   end subroutine Data_Book_Clear!}}}

   subroutine Gen_Array_SetCapacity(L, C)!{{{
   ! Set size of array
    Type (Gen_Array) :: L
    Integer :: C
    Double precision, dimension(:), pointer :: TmpValues

    if (L%length > 0) then
      if (C < L%length) stop 'Gen_Array_SetCapacity: smaller than length'
      allocate(TmpValues(L%length))
      TmpValues = L%Values(1:L%length)
      deallocate(L%Values)
      allocate(L%Values(C))
      L%Values(1:L%length) = TmpValues
      deallocate(TmpValues)
    else
      allocate(L%Values(C))
    end if
    L%Capacity = C
    
   end subroutine Gen_Array_SetCapacity!}}}

   subroutine Data_Set_SetCapacity(L, C)!{{{
    Type (Data_Set) :: L
    integer C
    double precision, dimension(:), pointer :: TmpItems_indep
    double precision, dimension(:,:), pointer :: TmpItems_dep
   
    if (L%n_indep > 0) then
      ! Already have things inside
      if (C < L%n_indep) stop 'Data_Set_SetCapacity: smaller than n_indep'
      allocate(TmpItems_indep(L%n_indep))
      allocate(TmpItems_dep(L%n_indep,L%n_dep))

      TmpItems_indep = L%Indep_Values(1:L%n_indep)
      TmpItems_dep = L%Dep_Values(1:L%n_indep,1:L%n_dep)

      deallocate(L%Indep_Values)
      deallocate(L%Dep_Values)

      allocate(L%Indep_Values(C))
      allocate(L%Dep_Values(C,L%n_dep))

      L%Indep_Values(1:L%n_indep) = TmpItems_indep
      L%Dep_Values(1:L%n_indep,1:L%n_dep) = TmpItems_dep

      deallocate(TmpItems_indep)
      deallocate(TmpItems_dep)
    else
      allocate(L%Indep_Values(C))
      allocate(L%Dep_values(C,L%n_dep))
    end if  
    L%Capacity = C
  
   end subroutine Data_Set_SetCapacity!}}}

   subroutine Data_Book_SetCapacity(L, C)!{{{
    Type (Data_Book) :: L
    integer C
    type(Data_Set), dimension(:), pointer :: TmpDatasets
   
    if (L%Count > 0) then
      ! Already have things inside
      if (C < L%Count) stop 'Data_Book_SetCapacity: smaller than Count'
      allocate(TmpDatasets(L%Count))

      TmpDatasets = L%Datasets(1:L%Count)

      deallocate(L%Datasets)
      allocate(L%Datasets(C))

      L%Datasets(1:L%Count) = TmpDatasets

      deallocate(TmpDatasets)
    else
      allocate(L%Datasets(C))
    end if  
    L%Capacity = C
  
   end subroutine Data_Book_SetCapacity!}}}

  subroutine Gen_Array_Add_Element(L, AVal)!{{{
    type (Gen_Array) :: L
    Double precision :: AVal

    if (L%length > 0) then 
      ! Array already exists, merely adding another element
      if (L%length == L%Capacity) call Gen_Array_SetCapacity(L, L%Capacity + Tab_blocksize)
    else ! nothing is present so start an array
      call Gen_Array_SetCapacity(L, Tab_blocksize) ! start off with 128 elements
    end if
    L%length = L%length + 1
    L%Values(L%length) = AVal
  end subroutine Gen_Array_Add_Element!}}}

   subroutine Data_Set_Add_Row(L, AVal, ARow)!{{{
    Type (Data_Set) :: L
    Double precision :: AVal
    Type (Gen_Array), intent(in) :: ARow
    integer :: rowid, i

    if (L%n_indep > 0) then
      ! Array already exists, merely adding another row
      if (L%n_dep > ARow%length) then
        Print *,'Dataset ',L%Name,' needs atleast ',L%n_dep,' elements.'
        stop
      end if
      if (L%n_indep == L%Capacity) call Data_Set_SetCapacity(L, L%Capacity + Tab_blocksize)
    else ! nothing is present so start an array
      L%n_dep = ARow%Length ! number of columns
      call Data_Set_SetCapacity(L, Tab_blocksize) ! Start off with 128 rows
    end if
    rowid = L%n_indep + 1 ! Adding one more to the values
    L%n_indep = rowid
    L%Indep_Values(rowid) = AVal
    Do i=1, L%n_dep
      L%Dep_Values(rowid,i) = ARow%Values(i)
    End do

   end subroutine Data_Set_Add_Row!}}}

  subroutine Data_Set_Add_Rowline(L, AInLine)!{{{
  ! Add a row to a table in an existing dataset from a string 
    Type(Data_Set) :: L
    character (LEN=*), intent(IN) :: AInLine

    integer :: Hashpos, lastpos, spos
    character (LEN=len(AInLine)) :: InLine, S
    Double precision AVal, AElement
    type(Gen_Array) :: ARow

    InLine=trim(adjustl(AInLine)) ! move all leading spaces to the end
    if (InLine=='' .or. InLine(1:1)=='#' .or. InLine(1:7)=='COMMENT' ) return

    Hashpos=scan(InLine,'#')
    if (Hashpos /= 0) then
      InLine  = InLine(1:Hashpos-1)
    end if

    lastpos=len_trim(InLine)
    if (lastpos>1) then
      if (InLine(1:1)=='''' .and. InLine(lastpos:lastpos)=='''') then
        InLine = InLine(2:lastpos-1)
      end if
    end if

    ! Read first element
    spos = Scan(InLine,' ')

    if (spos==0) then
      print *,'I need atleast two columns in dataset ',L%Name
      stop
    else
      S = trim(adjustl(InLine(1:spos-1)))
      InLine = trim(adjustl(InLine(spos+1:)))
      Read(S,*,err=40) AVal
    end if

    ! Define Gen_Array that contains values
    call Gen_Array_Init(ARow)

    ! Now scan through line and pass space separated values
    Do
      spos = Scan(InLine,' ')
      If (spos>1) then
        S = trim(adjustl(InLine(1:spos-1)))
        InLine = trim(adjustl(InLine(spos+1:)))
        Read(S,*,err=40) AElement
        call Gen_Array_Add_Element(ARow, AElement)
      Else if (InLine /= '') then
        S = trim(adjustl(InLine))
        Read(S,*,err=40) AElement
        call Gen_Array_Add_Element(ARow, AElement)
        exit
      Else
        exit
      End if
    End do

    ! now read the row into the dataset
    call Data_Set_Add_Row(L, AVal, ARow)
    ! Clear the Gen_Array created and free up memory
    call Gen_Array_Clear(ARow)
    return

 40 Print *,'Unable to read into double: ',S
    stop
  end subroutine Data_Set_Add_Rowline!}}}

  subroutine Data_Set_Add_Auxline(L, AInline)!{{{
  ! Adds aux array found on comment line to the dataset
    type(Data_Set) :: L
    character (LEN=*), intent(IN) :: AInLine

    integer :: lastpos, spos
    character (LEN=len(AInLine)) :: InLine, S
    Double precision AElement

    InLine=trim(adjustl(AInLine)) ! move all leading spaces to the end
    if (InLine=='') return
    lastpos=len_trim(InLine)
    if (lastpos>1) then
      if (InLine(1:1)=='''' .and. InLine(lastpos:lastpos)=='''') then
        InLine = InLine(2:lastpos-1)
      end if
    end if

    ! find where array of numbers starts
    Do
      spos = Scan(InLine,'-+0123456789.')
      if (spos==0) return
      ! Check if we found a sign, if not go ahead
      If ((InLine(spos:spos)/='-').and.(InLine(spos:spos)/='+')) exit
      If (len_trim(InLine)==spos) return ! no point moving ahead
      ! if we did find a sign, check that it precedes a number
      If (Scan(InLine(spos+1:),'0123456789')==1) exit
      ! it's some other sign, move ahead
      InLine = InLine(spos+1:)
    End do
    InLine = InLine(spos:)

    ! clear out existing aux, if any
    call Gen_Array_Clear(L%Aux_Values)

    ! Now scan through line and pass space separated values
    Do
      spos = Scan(InLine,' ') 
      If (spos>1) then
        S = trim(adjustl(InLine(1:spos-1)))
        InLine = trim(adjustl(InLine(spos+1:)))
        Read(S,*,err=40) AElement
        call Gen_Array_Add_Element(L%Aux_Values, AElement)
      Else if (InLine /= '') then
        S = trim(adjustl(InLine))
        Read(S,*,err=40) AElement
        call Gen_Array_Add_Element(L%Aux_Values, AElement)
        exit
      Else
        exit
      End if
    End do
    return

 40 Print *,'Unable to read into double: ',S
    stop
  end subroutine Data_Set_Add_Auxline!}}}

  function Data_Book_HasKey(Key) result(AValue)!{{{
   character (LEN=*), intent(IN) :: Key
   logical AValue

   AValue = Data_Book_HasKey_File(Def_Data, Key)

  end function Data_Book_HasKey!}}}

  function Data_Book_HasKey_File(L, Key) result(AValue)!{{{
    type(Data_Book), intent(in) :: L
    character (LEN=*), intent(IN) :: Key
    logical :: AValue
    integer :: i

    do i=1, L%Count
      if (L%Datasets(i)%Name == Key) then
        AValue = .true.
        return
      end if
    end do
    AValue = .false.     
    
  end function Data_Book_HasKey_File!}}}

   function Data_Book_IndexOf(AName) result (AValue)!{{{
     character(LEN=*), intent(in) :: AName
     integer :: AValue
     AValue = Data_Book_IndexOf_File(Def_Data, AName)
   end function Data_Book_IndexOf!}}}

   function Data_Book_IndexOf_File(L, Key) result (AValue)!{{{
     type(Data_Book) :: L
     character(LEN=*), intent(in) :: Key
     integer :: i, AValue

     do i=1, L%Count
       if (L%Datasets(i)%Name == Key) then
          AValue = i
          return
       end if
     end do
     AValue = -1
   end function Data_Book_IndexOf_File!}}}

   function Data_Book_NIndep(AName) result (AValue)!{{{
     character(LEN=*), intent(in) :: AName
     integer :: AValue
     AValue = Data_Book_NIndep_File(Def_Data, AName)
   end function Data_Book_NIndep!}}}

   function Data_Book_NIndep_File(L, Key) result (AValue)!{{{
     type(Data_Book) :: L
     character(LEN=*), intent(in) :: Key
     integer :: i, AValue

     if (.not. Data_Book_HasKey_File(L,Key)) then
       print *,"No data by the name of ",Key
       stop
     end if

     i = Data_Book_IndexOf_File(L, Key)
     AValue = L%Datasets(i)%n_indep
   end function Data_Book_NIndep_File!}}}

   function Data_Book_NDep(AName) result (AValue)!{{{
     character(LEN=*), intent(in) :: AName
     integer :: AValue
     AValue = Data_Book_NDep_File(Def_Data, AName)
   end function Data_Book_NDep!}}}

   function Data_Book_NDep_File(L, Key) result (AValue)!{{{
     type(Data_Book) :: L
     character(LEN=*), intent(in) :: Key
     integer :: i, AValue

     if (.not. Data_Book_HasKey_File(L,Key)) then
       print *,"No data by the name of ",Key
       stop
     end if

     i = Data_Book_IndexOf_File(L, Key)
     AValue = L%Datasets(i)%n_dep
   end function Data_Book_NDep_File!}}}

  subroutine Data_Open(filename, unit_id,  AName, Aux_line)!{{{
     character (LEN=*), intent(IN) :: filename
     integer, intent(IN) :: unit_id
     character (LEN=*), intent(IN) :: AName
     integer, optional, intent(IN) :: Aux_line
     logical aerror
          
     if (present(Aux_line)) then
      call Data_Open_File(Def_Data,filename,unit_id,AName,aerror,Aux_line)
     else
      call Data_Open_File(Def_Data,filename,unit_id,AName,aerror)
     end if

     if (aerror) then
        write (*,*) 'Data_Open: Error opening file ' // trim(filename)
        stop
     end if
  end subroutine Data_Open!}}}

  recursive subroutine Data_Open_File(L, filename, unit_id,             &!{{{
     &                               AName, error, Aux_line)
    Type(Data_Book) :: L
    character (LEN=*), intent(IN) :: filename
    integer, intent(IN) :: unit_id
    character (LEN=*), intent(IN) :: AName
    logical, intent(OUT) :: error
    integer, optional, intent(IN) :: Aux_line

    character (LEN=TAB_MAX_ROW_LEN) :: InLine
    integer i

    if (Data_Book_HasKey_File(L, AName)) then
      print *,"Ive already read in data for ",Aname
      return
    end if

    if (L%Count > 0) then
      ! There's already other data in there
      if (L%Count == L%Capacity) call Data_Book_SetCapacity(L, L%Capacity + Tab_blocksize)
    else ! Nothing is present so start 
      call Data_Book_SetCapacity(L, Tab_blocksize) ! start off with 128 files
    end if
    L%Count = L%Count + 1
    ! Adding into L%Datasets(L%Count)
    call Data_Set_Init(L%Datasets(L%Count))
    L%Datasets(L%Count)%Name = AName
             
    open(unit=unit_id,file=filename,form='formatted',status='old', err=500)

    i = 0
    do 
      read (unit_id,'(a)',end=400) InLine
      i = i + 1
      if (InLine == 'END') exit;
      if (InLine /= '') then
       call Data_Set_Add_Rowline(L%Datasets(L%Count),InLine)
       if (present(Aux_line)) then
         if (i==Aux_line) then
           call Data_Set_Add_Auxline(L%Datasets(L%Count),InLine)
         end if
       end if
      end if
    end do

400 close(unit_id)
    error=.false.
    If (IO_VERBOSE) Print *,"Read in ",L%Datasets(L%Count)%n_indep,     &
     &                      " independent parameters and ",             &
     &                      L%Datasets(L%Count)%n_dep,                  &
     &                      " dependent parameters from ",trim(filename) 
    return

500 error=.true.
    return
  end subroutine Data_Open_File!}}}

  subroutine Data_Read( Key,Indep,Indep_in_size,Indep_out_size,Dep,     &!{{{
     & Dep_in_i_size, Dep_in_j_size, Dep_out_j_size,                    &
     & Aux,Aux_in_size,Aux_out_size )
     ! This subroutine accepts a key, and reads out the data labeled 
     ! by the key. The independent variable is finally in Indep, and the 
     ! dependent variables' data is in Dep. If Aux is asked for, it is 
     ! placed in Aux. The in sizes are for bound checking, and the out
     ! sizes hold the data size

     character (LEN=*), intent(IN) :: Key
     double precision, dimension(:), intent(OUT) :: Indep
     double precision, dimension(:,:), intent(OUT) :: Dep
     integer, intent(IN) :: Indep_in_size, Dep_in_i_size, Dep_in_j_size
     integer, intent(OUT) :: Indep_out_size, Dep_out_j_size

     double precision, dimension(:), optional, intent(OUT) :: Aux
     integer, optional, intent(IN) :: Aux_in_size
     integer, optional, intent(OUT) :: Aux_out_size

     if (present(Aux) .and. present(Aux_in_size) .and. present(Aux_out_size)) then
     call Data_Read_File(Def_Data,Key,Indep,Indep_in_size,              &
     & Indep_out_size,Dep,Dep_in_i_size, Dep_in_j_size,                 &
     & Dep_out_j_size,Aux,Aux_in_size,Aux_out_size )
     else
     call Data_Read_File(Def_Data,Key,Indep,Indep_in_size,              &
     & Indep_out_size,Dep,Dep_in_i_size, Dep_in_j_size,                 &
     & Dep_out_j_size)
     end if

   end subroutine Data_Read!}}}

  subroutine Data_Read_File(L,Key,Indep,Indep_in_size,Indep_out_size,   &!{{{
     & Dep,Dep_in_i_size, Dep_in_j_size, Dep_out_j_size,                &
     & Aux,Aux_in_size,Aux_out_size )

     type(Data_Book), intent(IN) :: L
     character (LEN=*), intent(IN) :: Key
     double precision, dimension(:), intent(OUT) :: Indep
     double precision, dimension(:,:), intent(OUT) :: Dep
     integer, intent(IN) :: Indep_in_size, Dep_in_i_size, Dep_in_j_size
     integer, intent(OUT) :: Indep_out_size, Dep_out_j_size

     double precision, dimension(:), optional, intent(OUT) :: Aux
     integer, optional, intent(IN) :: Aux_in_size
     integer, optional, intent(OUT) :: Aux_out_size

     integer :: m, i, j

     if (.not. Data_Book_HasKey_File(L,Key)) then
       print *,"No data by the name of ",Key
       stop
     end if

     m = Data_Book_IndexOf_File(L, Key)

     if ((L%Datasets(m)%n_indep > Indep_in_size) .or.                   &
     &   (L%Datasets(m)%n_indep > Dep_in_i_size) .or.                   &
     &   (L%Datasets(m)%n_dep > Dep_in_j_size)) then
       print *,"input array not big enough to hold data for ",Key
       print *,"increase max_data_size"
       stop
     end if

     Indep_out_size = L%Datasets(m)%n_indep
     Dep_out_j_size = L%Datasets(m)%n_dep

     Do i=1, Indep_out_size
       Indep(i) = L%Datasets(m)%Indep_Values(i)
     End do

     Do j=1, Dep_out_j_size
       Do i=1, Indep_out_size
         Dep(i,j) = L%Datasets(m)%Dep_Values(i,j)
       End do
     End do

     if (present(Aux) .and. present(Aux_in_size) .and. present(Aux_out_size)) then
       if (L%Datasets(m)%Aux_Values%length > Aux_in_size) then
         print *,"input aux array not big enough to hold data for ",Key
         print *,"increase max_data_size"
         stop
       end if
       Aux_out_size = L%Datasets(m)%Aux_Values%length
       Do i=Aux_out_size, 1, -1
         Aux(i) = L%Datasets(m)%Aux_Values%Values(i)
       End do
     end if
     return

   end subroutine Data_Read_File!}}}

end module Data_Repo!}}}
!---------------------------------------------------------------------
