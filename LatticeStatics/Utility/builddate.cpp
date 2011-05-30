const char* builddate();

const char* builddate()
{
#ifdef BUILD_DATE
   return BUILD_DATE;
#else
   return "Unknown";
#endif
}
