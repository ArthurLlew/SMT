#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
    printf("Hello, World!\n");
    for (int i = 0; i < argc; i++)
    {
        int a = atoi(argv[i]);
        printf("%d\n",a);
    }
    return 0;
}