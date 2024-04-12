
//thiet lap module DS1307
//DFRobot_DS1307 DS1307;
LiquidCrystal_I2C lcd(0x27,16,2);  // set the LCD address to 0x27 for a 16 chars and 2 line display
unsigned long prvReadTime = 0;        // will store last time LED was updated

// constants won't change:
const long iReadTime = 1000;           // thoi gian doc DS1307 la 1s

tmElements_t tm;
#define mucnuoc_PIN A0                      //chan doc muc nuoc
int mucnuoc = 0;
int mucnuocthap =30;
int mucnuoccao = 50;
int gio,phut;
unsigned long prvHienThiLCD = 0;        // will store last time LED was updated
const long iHienThiLCD = 100;           // thoi gian hien thi LCD

unsigned long prvBatMayBom = 0;        // will store last time LED was updated
const long iBatMayBom = 100;           // thoi gian hien thi LCD
#define maybom_PIN  9                  //dinh nghia chan may bom
#define led_PIN     10                 //dinh nghia chan Led chieu sang

boolean trangThaiMayBom =false;
boolean trangThaiLed = false;
boolean b_thucHienMotLan=false;
#define MODE_PIN    4
Bounce btnmode = Bounce();
#define SET_PIN     3
Bounce btnset = Bounce();
#define TANG_PIN    6
Bounce btntang = Bounce();
#define GIAM_PIN    5
Bounce btngiam = Bounce();
#define TEST_PIN 2
Bounce btntest = Bounce();
#define HUY_PIN     7
Bounce btnhuy = Bounce();
