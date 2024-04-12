#include <Bounce2.h>

#include <Wire.h> 
#include <LiquidCrystal_I2C.h>
//#include <DFRobot_DS1307.h>
#include <TimeLib.h>
#include <DS1307RTC.h>
#include "thamso.h"
#include "program.h"


void setup() {
    lcd.init();                      // initialize the lcd 
    lcd.init();
    lcd.backlight();                //sang den nen LCD
    lcd.setCursor(0,0);
    lcd.print("CDIO5-MN:");
    Serial.begin(115200);           //mo cong Serial
    pinMode(maybom_PIN,OUTPUT);
   // digitalWrite(maybom_PIN,LOW);        //thiet lap ngo ra
     pinMode(led_PIN,OUTPUT);
     //digitalWrite(led_PIN,LOW);             //thiet lap ngo ra
     
     
     digitalWrite(led_PIN,LOW); 
     digitalWrite(maybom_PIN,LOW);
      btnmode.attach(MODE_PIN,INPUT_PULLUP); // Attach the debouncer to a pin with INPUT_PULLUP mode
  btnmode.interval(25);
  btnset.attach(SET_PIN,INPUT_PULLUP); // Attach the debouncer to a pin with INPUT_PULLUP mode
  btnset.interval(25);
  btntang.attach(TANG_PIN,INPUT_PULLUP); // Attach the debouncer to a pin with INPUT_PULLUP mode
  btntang.interval(25);
  btngiam.attach(GIAM_PIN,INPUT_PULLUP); // Attach the debouncer to a pin with INPUT_PULLUP mode
  btngiam.interval(25);
  btnhuy.attach(HUY_PIN,INPUT_PULLUP); // Attach the debouncer to a pin with INPUT_PULLUP mode
  btnhuy.interval(25);
  btnset.attach(SET_PIN,INPUT_PULLUP); // Attach the debouncer to a pin with INPUT_PULLUP mode
  btnset.interval(25);
  
   delay(1000);
    
    //-----ham de setup thoi gian ban dau
    /*
     while( !(DS1307.begin()) ){
        Serial.println("Communication with device failed, please check connection");
        delay(3000);
    }
    Serial.println("Begin ok!");
    DS1307.setTypeTime(DS1307.eYR, 2023);
    uint16_t setTimeBuff[7] = {5, 15, 14, 6, 21, 6, 2023};
    DS1307.setTime(setTimeBuff);
     DS1307.start();
    */
}

void loop() {
    readtime();         //HAM DOC THOI GIAN TU DS1307
    hienthiLCD();
    hambatmaybom();
    hamChieuSang(2);
    readbutton();
//    hamkiemtraden();
   // hamkiemtramaybom();
    //testmaybom();
   // docMucNuoc();
}
