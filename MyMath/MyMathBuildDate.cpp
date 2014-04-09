char* MyMathBuildDate();

char* MyMathBuildDate()
{
#ifdef BUILD_DATE
  return (char*) BUILD_DATE;
#else
  return (char*) "Unknown";
#endif
}
