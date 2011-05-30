char* MyMathBuildDate();

char* MyMathBuildDate()
{
#ifdef BUILD_DATE
   return BUILD_DATE;
#else
   return "Unknown";
#endif
}
