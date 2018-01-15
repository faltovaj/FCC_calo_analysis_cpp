#include "Decoder.h"
#include <boost/algorithm/string.hpp>
#include <boost/lexical_cast.hpp>
#include <vector>
#include <iostream>

Decoder::Decoder(const std::string& aDescr): m_description(aDescr), m_currentValue(0) {

  // create a vector of fields:
  std::vector<std::string> fields;
  boost::split(fields, aDescr, boost::is_any_of(","));
  unsigned long long offset = 0;
  for(auto s: fields) {
    std::cout<<s<<"\t"<<std::endl;
    std::vector<std::string> onefield;
    boost::split(onefield, s, boost::is_any_of(":"));
    uint w =  abs(boost::lexical_cast<int>(onefield[1]));
    Field field;
    field.name = onefield[0];
    field.width = w;
    field.offset = offset;
    field.isSigned = (boost::lexical_cast<int>(onefield[1]) < 0);
    field.mask = ( ( 0x0001LL << field.width ) - 1 ) << field.offset;
    if( field.isSigned ){
      field.minVal =  - ( ( 1LL << ( field.width - 1 ) ) - 1);
      field.maxVal =  ( 1LL << ( field.width - 1 ) );
    } else {
      field.minVal = 0;
      field.maxVal = ( 1LL << field.width ) - 1;
    }
    m_fields.insert(make_pair(field.name, field));
    std::cout << "\tfield name:" << field.name
              << "\tfield width:" << field.width
              << "\tfield signed:" << field.isSigned
              << "\tfield offset:" << std::hex <<field.offset << std::dec
              << "\tfield min value:" << field.minVal
              << "\tfield max value:" << field.maxVal
              << "\tfield mask:" << std::hex << field.mask << std::dec
              << std::endl;
    offset += w;
  }

}

unsigned long long Decoder::operator[](const std::string& aName) const {
  return value(aName);
}


unsigned long long Decoder::value(const std::string& aName) const {
  // Not checking if field is signed
  // TODO
  return (m_currentValue & m_fields.find(aName)->second.mask) >> m_fields.find(aName)->second.offset;
}

unsigned long long Decoder::value(const std::string& aName, long long aValue) const {
  // Not checking if field is signed
  // TODO
  return (aValue & m_fields.find(aName)->second.mask) >> m_fields.find(aName)->second.offset;
}

void Decoder::setValue(long long aValue) {
  m_currentValue = aValue;
}

double  Decoder::segmentationPosition(long long aCellValue, double aCellSize, double aCellOffset) const {
  return (aCellValue * aCellSize + aCellOffset);
}
