#include <stdio.h>
#include <stdio.h>

int main (int argc, char** argv)
{
  FILE *f1 = fopen(argv[1], "r");
  FILE *f2 = fopen(argv[2], "r");
  int cnt1 = 0, cnt2 = 0;
  int c1, c2; 
 
  while((c1=fgetc(f1))!=EOF){
    cnt1++;
    c2=fgetc(f2);
    if(c1 != c2){
      printf("char pos %i, char 1 %c, char 2 %c\n", cnt1,c1,c2);
}
}
}
  
