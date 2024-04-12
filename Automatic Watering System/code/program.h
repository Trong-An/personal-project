#include <Bounce2.h>

void hamChieuSang(int _thoigian);

void print2digits(int number) {
  if (number >= 0 && number < 10) {
    Serial.write('0');
  }
  Serial.print(number);
}
void readtime()         //doc thoi gian tu DS1307
{
    unsigned long cReadTime = millis();

  if (cReadTime - prvReadTime >= iReadTime) 
  {
        prvReadTime = cReadTime;
        

    if (RTC.read(tm)) {
        Serial.print("Ok, Time = ");
        print2digits(tm.Hour);
        Serial.write(':');
        print2digits(tm.Minute);
        Serial.write(':');
        print2digits(tm.Second);
        Serial.print(", Date (D/M/Y) = ");
        Serial.print(tm.Day);
        Serial.write('/');
        Serial.print(tm.Month);
        Serial.write('/');
        Serial.print(tmYearToCalendar(tm.Year));
        Serial.println();

        
    } else 
    {
        if (RTC.chipPresent()) {
        Serial.println("The DS1307 is stopped.  Please run the SetTime");
        Serial.println("example to initialize the time and begin running.");
        Serial.println();
        } else {
        Serial.println("DS1307 read error!  Please check the circuitry.");
        Serial.println();
        }
        delay(9000);
    }
    
  }
}
void hienthiLCD()
{
   unsigned long cHienThiLCD = millis();

  if (cHienThiLCD - prvHienThiLCD >= iHienThiLCD) 
  {
    prvHienThiLCD = cHienThiLCD;
    
    //doc muc nuoc
    mucnuoc = map(analogRead(mucnuoc_PIN),0,1024,0,100);
    lcd.setCursor(0,0);
    char outputarr[16];
    sprintf(outputarr, "CDIO5-MN:%3d",mucnuoc);
    lcd.print(outputarr);
    
    //hien thi gio phut giay
    char outputarr2[16];
    sprintf(outputarr2, "%02d:%02d:%02d  %1d %1d",tm.Hour,tm.Minute,tm.Second,trangThaiMayBom,trangThaiLed);
    lcd.setCursor(0,1);
    lcd.print(outputarr2);

  }

}
/*void testden(){
  int buttonStatus = digitalRead(sw1_PIN);
  if(buttonStatus == HIGH){
    digitalWrite(led_PIN,HIGH);
    }
  }*/
/*void hamkiemtraden(){
    int buttonStatus = digitalRead(sw1_PIN);
  if(buttonStatus == LOW){
     b_thucHienMotLan = true;
     digitalWrite(led_PIN,HIGH);
    }
    else{
       b_thucHienMotLan = true;
           digitalWrite(led_PIN,LOW);
      }
    }
void hamkiemtramaybom(){
    int buttonStatus = digitalRead(sw2_PIN);
  if(buttonStatus == LOW){
     trangThaiMayBom = true;
     digitalWrite(maybom_PIN,HIGH); 
    }
    else{
       trangThaiMayBom = true;
       digitalWrite(maybom_PIN,LOW); 
      }
    }*/
void readbutton()
{
    //btnmode.update(); // Update the Bounce instance
   
   // if ( btnmode.fell() ) {  // Call code if button transitions from HIGH to LOW
       // bip(2,100);
   // }

     btntest.update();
    if(btntest.fell()){
      digitalWrite(led_PIN,HIGH);
      char outputarr[16];
      lcd.clear();
      sprintf(outputarr,"MOTOR ON");
      lcd.setCursor(0,1);
      }
    
    btnset.update(); // Update the Bounce instance
   
    if ( btnset.fell() ) {  // Call code if button transitions from HIGH to LOW
        //if(mode == modeThietLapNguong)
        //{
            int modesetnguong=0;
            while(1)
            {
                if(modesetnguong==0)
                {
                        char outputarr1[16];
                        sprintf(outputarr1, ">Nguong Thap:%03d",mucnuocthap);
                        lcd.setCursor(0,0);lcd.print(outputarr1);
                        char outputarr2[16];
                            sprintf(outputarr2, "Nguong Cao:  %03d",mucnuoccao);
                        lcd.setCursor(0,1);lcd.print(outputarr2);

                        btntang.update(); // Update the Bounce instance
        
                        if ( btntang.fell() ) 
                        {  
                            mucnuocthap++;
                        }
                        btngiam.update(); // Update the Bounce instance
   
                        if ( btngiam.fell() ) {  // Call code if button transitions from HIGH to LOW
                            mucnuocthap--;
                        }
                         
                         btnset.update(); // Update the Bounce instance
   
                        if ( btnset.fell() )
                        {
                            modesetnguong=1;
                        }
                        btnhuy.update(); // Update the Bounce instance
   
                        if ( btnhuy.fell() ) {  // Call code if button transitions from HIGH to LOW
                            break;
                        }
  
                }
                if(modesetnguong==1)
                {
                        char outputarr1[16];
                        sprintf(outputarr1, "Nguong Thap: %03d",mucnuocthap);
                        lcd.setCursor(0,0);lcd.print(outputarr1);
                        char outputarr2[16];
                            sprintf(outputarr2, ">Nguong Cao: %03d",mucnuoccao);
                        lcd.setCursor(0,1);lcd.print(outputarr2);

                        btntang.update(); // Update the Bounce instance
        
                        if ( btntang.fell() ) 
                        {  
                            mucnuoccao++;
                        }
                        btngiam.update(); // Update the Bounce instance
   
                        if ( btngiam.fell() ) {  // Call code if button transitions from HIGH to LOW
                            mucnuoccao--;
                        }
                         
                         btnset.update(); // Update the Bounce instance
   
                        if ( btnset.fell() )
                        {
                            modesetnguong=0;
                        }

                        btnhuy.update(); // Update the Bounce instance
   
                        if ( btnhuy.fell() ) {  // Call code if button transitions from HIGH to LOW
                            break;
                        }
                    
                }
                                    
            }
    }
    btnhuy.update();
    if ( btnhuy.fell() ) {  // Call code if button transitions from HIGH to LOW
       
       
           int modesetthoigian=0;
            while(1)
            {
               if(modesetthoigian==0)
                {
                        char outputarr1[16];
                        sprintf(outputarr1, ">Gio: %02d    ",gio);
                        lcd.setCursor(0,0);lcd.print(outputarr1);
                        char outputarr2[16];
                            sprintf(outputarr2, "Phut:  %02d   ",phut);
                        lcd.setCursor(0,1);lcd.print(outputarr2);

                        btntang.update(); // Update the Bounce instance
        
                        if ( btntang.fell() ) 
                        {  
                            gio++;
                        }
                        btngiam.update(); // Update the Bounce instance
   
                        if ( btngiam.fell() ) {  // Call code if button transitions from HIGH to LOW
                            gio--;
                        }
                         
                         btnset.update(); // Update the Bounce instance
   
                        if ( btnset.fell() )
                        {
                            modesetthoigian=1;
                        }
                        btnhuy.update(); // Update the Bounce instance
   
                        if ( btnhuy.fell() ) {  // Call code if button transitions from HIGH to LOW
                            break;
                        }
  
                }
                if(modesetthoigian==1)
                {
                        char outputarr1[16];
                        sprintf(outputarr1, "Gio: %02d   ",gio);
                        lcd.setCursor(0,0);lcd.print(outputarr1);
                        char outputarr2[16];
                            sprintf(outputarr2, ">Phut: %02d   ",phut);
                        lcd.setCursor(0,1);lcd.print(outputarr2);

                        btntang.update(); // Update the Bounce instance
        
                        if ( btntang.fell() ) 
                        {  
                            phut++;
                        }
                        btngiam.update(); // Update the Bounce instance
   
                        if ( btngiam.fell() ) {  // Call code if button transitions from HIGH to LOW
                            
                            phut--;
                        }
                         
                         btnset.update(); // Update the Bounce instance
   
                        if ( btnset.fell() )
                        {
                            
                            modesetthoigian=0;
                        }

                        btnhuy.update(); // Update the Bounce instance
   
                        if ( btnhuy.fell() ) {  // Call code if button transitions from HIGH to LOW
                           
                            break;
                        }
                }
               

            }
    }


        
}    
void hambatmaybom()
{
    if(mucnuoc<50)
    {
        if(!trangThaiMayBom)
        {
            trangThaiMayBom = true;
            digitalWrite(maybom_PIN,HIGH); 
            while(mucnuoc<=60)
            {
                readtime();
                //mucnuoc = map(analogRead(mucnuoc_PIN),0,1024,0,100);
                hienthiLCD();
                hamChieuSang(2);
            } 
             trangThaiMayBom = false;
            digitalWrite(maybom_PIN,LOW); 
        }
    }
}
void hamChieuSang(int _thoigian)
{
    if(tm.Minute%_thoigian==0 && !b_thucHienMotLan)          //sau 2 phut bat tat 1 lan
    {
           //trangThaiLed = !trangThaiLed;
           b_thucHienMotLan = true;
           digitalWrite(led_PIN,HIGH);
    }
    if(tm.Minute%_thoigian !=0)
    {
        b_thucHienMotLan = false;
         digitalWrite(led_PIN,LOW);
       // trangThaiLed=false;
    }
}
