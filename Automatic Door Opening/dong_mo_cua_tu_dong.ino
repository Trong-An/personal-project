
#include <Keypad.h>
#include <LiquidCrystal.h>
#include <EEPROM.h>
const int rs = 13, en = 12, d4 = 11, d5 = 10, d6 = 9, d7 = 8;
LiquidCrystal lcd(rs, en, d4, d5, d6, d7);
char chuoi[20];


// khao bao phim ma tran
const byte rows = 4; //so hang
const byte columns = 4; //so cot

int holdDelay = 200; //time tao tre nhan nut, tranh nhieu nut nhan
int n = 30; // thoi gian nhan giu nut
int state = 0; // state =0 ko nhan nut,state =1 nhan thoi gian nhan , state = 2 nhan giu nut
int key = 0; // gia tri nut khi chua dc nhan
char temp=0; // bien de doc gia tri nut 
//dinh nghia cac nut nhan tra ve
char keys[rows][columns] =
{
  {'7', '8', '9', '%'},
  {'4', '5', '6', 'x'},
  {'1', '2', '3', '-'},
  {'A', '0', '=', '+'},
};

 //Cach noi chan voi Arduino
 
byte rowPins[rows] = {0, 1, 2, 3}; 
byte columnPins[columns] = {4, 5, 6, 7};

//khai bao khoi tao keypad
Keypad keypad = Keypad(makeKeymap(keys), rowPins, columnPins, rows, columns);

  // put your setup code here, to run once:
#define da_dong 0
#define da_mo 1
#define o_giua 2
#define tu_dong 1
#define dk_tay 2

#define time_dong_cua 2000
#define time_mo_cua 2000

//const int coi_bao_dong =  A2;                          
//const int cam_bien = A4;  
const int dong_cua = A0; 
const int mo_cua =  A1;                         

int trang_thai_cb = 0;   
int trang_thai_coi = LOW; 
int trang_thai_cua=da_mo;
int thoi_gian = 0;
long delay_Start;
int k;

char che_do=0;
char cai_dat=0;
char lua_chon=0;
char pass_ok;
char check_pass;
char so_ki_tu,ki_tu;
unsigned char so_lan_con_lai;
char pass_temp[4];
char pass[4];
char pass_moi[4];

void Pir_Power_On();                             // cho module cam bien  khoi dong
void dong_cua_lai(int tg);    
void mo_cua_ra(int tg);  
void dung_motor();   
void nhap_pass_moi();
void xac_nhan_pass();
void hien_thi_lua_chon();    
void ctr_tu_dong();    
void ctr_dk_tay();              
void setup() {
  lcd.begin(16, 2);
  che_do = EEPROM.read(0x10);
  if(che_do>2 || che_do==0)
    che_do=dk_tay;
  for(so_ki_tu=0;so_ki_tu<4;so_ki_tu++){
    pass[so_ki_tu] = EEPROM.read(so_ki_tu);
    if(pass[so_ki_tu]<0||pass[so_ki_tu]>9){
      pass[so_ki_tu] =0;
    }
  }
 // EEPROM.write(0x11,5);
  so_lan_con_lai=EEPROM.read(0x11);
  if(so_lan_con_lai>5)
    so_lan_con_lai=5;
  else if(so_lan_con_lai==0)
    so_lan_con_lai=1;
  pinMode(coi_bao_dong, OUTPUT);
  pinMode(dong_cua, OUTPUT);
  pinMode(mo_cua, OUTPUT); 
  pinMode(cam_bien, INPUT);                   
//  Pir_Power_On();                               // cho module cam bien  khoi dong
 delay_Start = millis(); 
 thoi_gian=0;

 
}

void loop() {
  if(che_do==tu_dong){
    lcd.setCursor(0, 0);
    lcd.print(" CHE DO TU DONG");
    ctr_tu_dong();
  }
  else if(che_do==dk_tay){
    lcd.setCursor(0, 0);
    lcd.print(" CHE DO D.K TAY");
    ctr_dk_tay();
  }
  quet_phim();
  if(key==67) { // nhan giu nut menu
    hien_thi_lua_chon();
  }
  // put your main code here, to run repeatedly:
  
}

void bao_dong_coi(int lan_bao,int thoi_gian_bao){
  for(int i=0;i<lan_bao;i++)
  {
     digitalWrite(coi_bao_dong, HIGH);  
     delay(thoi_gian_bao);                               
     digitalWrite(coi_bao_dong, LOW);
     delay(thoi_gian_bao);  
  }
}
void Pir_Power_On(){
  delay(900);                               // cho modum pir on dinh dong xong
  //bao khoi dong xong
  bao_dong_coi(1,500);                  
}

void dong_cua_lai(){ 
  digitalWrite(mo_cua, LOW); 
  delay(1); 
  digitalWrite(dong_cua, HIGH); 
  trang_thai_cua=o_giua;                                                  
}
void mo_cua_ra(){
  digitalWrite(dong_cua, LOW); 
  delay(1); 
  digitalWrite(mo_cua, HIGH);  
   trang_thai_cua=o_giua;         
}
void dung_motor(){ 
  digitalWrite(dong_cua, LOW);  
  digitalWrite(mo_cua, LOW);                                                                     
}

void quet_phim(){ 
  temp = keypad.getKey();
  if ((int)keypad.getState() ==  PRESSED) {
    if (temp != 0) {
      key = temp;
    }
  } 

 if ((int)keypad.getState() ==  HOLD) {
    if (state==0){
        key = key+n;
        state=n;
    }
    delay(holdDelay);
  }

if ((int)keypad.getState() ==  RELEASED) {
    state = 0;// chuan bi cho lan nhan de nut tiep theo
    key=0;
    
  }
  delay(100);
                                                                     
}

void hien_thi_lua_chon(){ 
   lcd.clear();
    cai_dat=1; // set vao che do cai dat
    lua_chon=2;// lua chon cai dat che do, 1 la lua chon cai dat mat khau
    while(cai_dat==1){ // hien thi 1 so lua chon cho ng dung lua chon
      quet_phim();
      if(key==95) { // nhan giu nut on/ac
        lua_chon--;
        if(lua_chon==0)
          lua_chon=2;
        key=0;
      }
      else if(key==91) { // nhan giu nut =
        lua_chon++;
        if(lua_chon>2)
          lua_chon=1;
        key=0;
      }
      else if(key==150) { // nhan giu nut x => thoat
        lcd.clear();
        lua_chon=0;
        cai_dat=0;
        key=0;
      }
      else if(key==75) { // nhan giu nut ok 
        key=0;
        lcd.clear();
        lcd.setCursor(0, 0);
        lcd.print(" NHAP MAT KHAU! ");
        xac_nhan_pass();
        if(lua_chon==2){
          lcd.clear();
          lua_chon=4;// di toi cai dat che do
          cai_dat=2;// lua chon che do hoat dong
        }
        else if(lua_chon==1){
          lcd.clear();
          lua_chon=3;// di toi cai dat mk
          cai_dat=3;// thay doi mat khau 
        }
        
      }
 
      if(lua_chon==2){ // hien thi cai dat che do
        lcd.setCursor(0, 0);
        lcd.print("=> C.DAT CHE DO ");
        lcd.setCursor(0, 1);
        lcd.print("   C.DAT M.KHAU ");
      }
      else if(lua_chon==1){ // hien thi cai dat mat khau
        lcd.setCursor(0, 0);
        lcd.print("   C.DAT CHE DO ");
        lcd.setCursor(0, 1);
        lcd.print("=> C.DAT M.KHAU ");
      }
      
    while(cai_dat==2){ // bat che do tu dong hoac dk tay
      quet_phim();
      if(key==95) { // nhan giu nut on/ac
        lua_chon--;
        if(lua_chon==3)
          lua_chon=5;
        key=0;
      }
      else if(key==91) { // nhan giu nut =
        lua_chon++;
        if(lua_chon>5)
          lua_chon=4;
        key=0;
      }
      else if(key==150) { // nhan giu nut x => thoat
        lcd.clear();
        lua_chon=0;
        cai_dat=0;
        key=0;
        pass_ok=0;
      }
      else if(key==75) { // nhan giu nut ok 
        lcd.clear();
        if(lua_chon==4){
         che_do=tu_dong;
        }
        else if(lua_chon==5){
          che_do=dk_tay;
        }
        EEPROM.write(0x10,che_do);
        lua_chon=0;
        cai_dat=0;
        key=0;
        pass_ok=0;
      }
      
      if(lua_chon==4){ //  bat che do TU DONG
          lcd.setCursor(0, 0);
          lcd.print("=>CHE DO TU DONG");
          lcd.setCursor(0, 1);
          lcd.print("  CHE DO D.K TAY");
       }
       else if(lua_chon==5){ //   bat che do DK TAY
          lcd.setCursor(0, 0);
          lcd.print("  CHE DO TU DONG");
          lcd.setCursor(0, 1);
          lcd.print("=>CHE DO D.K TAY");
       }
        
    }

    while(cai_dat==3){ // thay doi mk 
        lcd.setCursor(0, 0);
        lcd.print(" NHAP PASS MOI! ");
        so_ki_tu=0;
        sprintf(pass_moi,"    ");// reset pass_temp
        nhap_pass_moi();
    }
  }                                                                   
}
void nhap_pass_moi(){

  char xac_nhan=0;
  char pass_moi_ok=0;
  so_ki_tu=0;
  while(so_ki_tu<4){
     quet_phim();
  
    if(key==150) { // nhan giu nut x => thoat
        lcd.clear();
        lua_chon=2;
        cai_dat=1;
        so_ki_tu=5;//thoat while
        key=0;
        pass_ok=0;
    }
    if( (key>47 && key<58) || (key>77 && key<88) ){ // du lieu nhap vao la so
         so_ki_tu++;
         lcd.setCursor(5+so_ki_tu, 1);
         lcd.print("*");
         if(xac_nhan==0){
            pass_moi[so_ki_tu-1]=key-48;
            if(key>57){
             pass_moi[so_ki_tu-1]=key-n-48;
            }
         }
         else if(xac_nhan==1){
            pass_temp[so_ki_tu-1]=key-48;
            if(key>57){
             pass_temp[so_ki_tu-1]=key-n-48;
            }
         }
         key=0;  
         if(so_ki_tu==4 && xac_nhan==0){
           lcd.setCursor(0, 0);
           lcd.print("XAC NHAN MK MOI"); 
           lcd.setCursor(0, 1);
           lcd.print("                "); 
           so_ki_tu=0;
           xac_nhan=1;
         }
         if(so_ki_tu==4 && xac_nhan==1){
            for(ki_tu=0;ki_tu<so_ki_tu;ki_tu++){
               if(pass_temp[ki_tu]==pass_moi[ki_tu]){ 
                  pass_moi_ok=1;
                  goto tiep;
               }
               else { // mk sai
                    sprintf(chuoi,"CON %d LAN NHAP",so_lan_con_lai);
                    lcd.setCursor(0, 1);
                    lcd.print("MK 2 LAN K.GIONG");
                    pass_moi_ok=0;
                    so_ki_tu=5;// thoat while;
                    delay(1000);
                    lcd.clear();
                    break;
               }
               tiep:;
            }
            if (pass_moi_ok==1){ // thuc hien het vong lap for ma k phat hien pass sai=> pass dung
                      lcd.setCursor(0, 1);
                      lcd.print("DA DOI MAT KHAU");
                      for(so_ki_tu=0;so_ki_tu<4;so_ki_tu++){
                         EEPROM.write(so_ki_tu,pass_moi[so_ki_tu]);
                         pass[so_ki_tu]=pass_moi[so_ki_tu];
                      }
                      so_ki_tu=5;// thoat vong lap while
            }
            
         }
    }
    if( key==43 ){ // nut +=> clr
         if(so_ki_tu>0){
            lcd.setCursor(5+so_ki_tu, 1);
            lcd.print(" ");
            if(xac_nhan==1)
              pass_temp[so_ki_tu]=10;
             else
              pass_moi[so_ki_tu]=10;
            so_ki_tu--;
         }
         key=0;
      }
  }
     
}
void xac_nhan_pass(){
    
    pass_ok=0;
    check_pass=0;
    so_ki_tu=0;
    sprintf(pass_temp,"    ");// reset pass_temp
    while(!check_pass){
      quet_phim();
      if( (key>47 && key<58) || (key>77 && key<88) ){ // du lieu nhap vao la so
         so_ki_tu++;
         lcd.setCursor(5+so_ki_tu, 1);
         lcd.print("*");
         pass_temp[so_ki_tu-1]=key-48;
         if(key>57){
           pass_temp[so_ki_tu-1]=key-n-48;
         }
         key=0;      
         if(so_ki_tu==4){
         
            for(ki_tu=0;ki_tu<so_ki_tu;ki_tu++){
               if(pass_temp[ki_tu]==pass[ki_tu]){ 
                  pass_ok=1;
                  goto next;
               }
               else { // mk sai
                    if(so_lan_con_lai>0){
                    so_lan_con_lai--;
                    sprintf(chuoi,"CON %d LAN NHAP",so_lan_con_lai);
                    lcd.setCursor(0, 1);
                    lcd.print(chuoi);
                    delay(1000);
                    lcd.setCursor(0, 1);
                    lcd.print("                ");
                  }
                  EEPROM.write(0x11,so_lan_con_lai);
                  if(so_lan_con_lai==0){ // nhap sai mk 5 lan
                      lcd.setCursor(0, 1);
                      lcd.print("XIN L.HE NHA SX!");
                  }
                  pass_ok=0;
                  so_ki_tu=0;
                  break;
               }
               next:;
            }
           if (pass_ok==1){ // thuc hien het vong lap for ma k phat hien pass sai=> pass dung
                      so_lan_con_lai=5;
                      check_pass=1;
                      so_ki_tu=0;
                      EEPROM.write(0x11,so_lan_con_lai);
            }
         }
      }
      
      
      if( key==43 ){ // nut +=> clr
         if(so_ki_tu>0){
            lcd.setCursor(5+so_ki_tu, 1);
            lcd.print(" ");
            pass_temp[so_ki_tu]=10;
            so_ki_tu--;
         }
         key=0;
      }

      if(key==150) { // nhan giu nut x => thoat
        lcd.clear();
        lua_chon=0;
        check_pass=1;
        cai_dat=0;
        key=0;
      }
      
    }
    
}

void ctr_tu_dong(){
  
  trang_thai_cb = digitalRead(cam_bien);
  if (trang_thai_cb == LOW && trang_thai_cua==da_mo) {  // neu khong co nguoi va cua da mo ra thi cho 1 khoang thoi gian de dong cua lai
       if((millis()- delay_Start) > 999){ 
          thoi_gian = thoi_gian+1;  
          delay_Start = millis(); // set start time
      }
      if(thoi_gian>5){ // da toi thoi gian can dong cua
          if(trang_thai_cua==da_mo){ // neu cua da mo  
              bao_dong_coi(2,250);   // thong bao chuan bi dong cua      
              lcd.setCursor(0, 1);
              lcd.print(" DANG DONG CUA! ");//  hien thi thong bao dang dong cua len man hinh lcd
              dong_cua_lai();
              // trong qua trinh dong cua ma co ng vao phai lap tuc mo cua ra ngay
              for(k=0;k<time_mo_cua;k++){
                delay(1);
                trang_thai_cb = digitalRead(cam_bien);
                if (trang_thai_cb == HIGH){ // luon luon kiem tra cam bien trong khi dong cua
                        dung_motor();
                        mo_cua_ra();
                        lcd.setCursor(0, 1);
                        lcd.print("  DANG MO CUA!  ");//  hien thi thong bao dang mo cua len man hinh lcd
                        delay(time_mo_cua+1000-k);// trc do da mo cua dc k ms=> khi mo cua chi can delay time_mo_cua-k => +1000 dam bao cua dong het
                        dung_motor();
                        trang_thai_cua=da_mo;
                        bao_dong_coi(4,50);// keu 4 tieng ngan là da mo
                        lcd.setCursor(0, 1);
                        lcd.print("   DA MO CUA!   ");//  hien thi thong bao da mo cua len man hinh lcd
                        thoi_gian=0;
                        delay_Start = millis(); // set start time
                        quet_phim();
                        if(key==67 ||key==37) { // nhan hoac nhan giu nut menu
                            hien_thi_lua_chon();
                        }  
                        break;
                }
              }
              if(k==time_mo_cua){
                trang_thai_cua=da_dong;
                dung_motor();
                bao_dong_coi(2,50);// keu 2 tieng ngan là da dong
                lcd.setCursor(0, 1);
                lcd.print("  DA DONG CUA!  ");//  hien thi thong bao da dong cua len man hinh lcd
                quet_phim();
                if(key==67 ||key==37) { // nhan hoac nhan giu nut menu
                  hien_thi_lua_chon();
                }
              }
             thoi_gian=0;
          } 
     }     
  }
  trang_thai_cb = digitalRead(cam_bien);
  if (trang_thai_cb == HIGH ) { // neu co nguoi di vao
       thoi_gian=0;
       delay_Start = millis(); // set start time
       if (trang_thai_cua!=da_mo){ // neu cua chua mo thi moi thuc hien 
             dung_motor();
             bao_dong_coi(1,250);  // thong bao chuan bi mo cua
             lcd.setCursor(0, 1);
             lcd.print("  DANG MO CUA!  ");//  hien thi thong bao dang mo cua len man hinh lcd 
             mo_cua_ra();
             delay(time_mo_cua);
             dung_motor();
             trang_thai_cua=da_mo;
             bao_dong_coi(4,50);// keu 4 tieng là da mo
             lcd.setCursor(0, 1);
             lcd.print("   DA MO CUA!   ");//  hien thi thong bao da mo cua len man hinh lcd
             thoi_gian=0;
             delay_Start = millis(); // set start time  
             quet_phim();
             if(key==67 ||key==37) { // nhan hoac nhan giu nut menu
                hien_thi_lua_chon();
             }      
       }
  } 
}

void ctr_dk_tay(){
  quet_phim();
   if(key==67 ||key==37) { // nhan hoac nhan giu nut menu
      hien_thi_lua_chon();
   }  
   else if(key==95) { // nhan giu nut on/ac => mo cua
      lcd.setCursor(0, 1);
      lcd.print("  DANG MO CUA!  ");//  hien thi thong bao dang mo cua len man hinh lcd
      mo_cua_ra();
   } 
   else if(key==91) { // nhan giu nut = => dong cua cua
      lcd.setCursor(0, 1);
      lcd.print(" DANG DONG CUA! ");//  hien thi thong bao dang dong cua len man hinh lcd
      dong_cua_lai();
   } 
   else if(key==0) { // nhan giu nut = => dong cua cua
      lcd.setCursor(0, 1);
      lcd.print("  DA DUNG CUA!  ");//  hien thi thong bao dang dong cua len man hinh lcd
      dung_motor();
   } 
}
