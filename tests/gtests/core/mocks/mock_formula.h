class Mock_Formula : public _Formula {
 public:
  MOCK_METHOD0(Initialize,
      void(void));
  MOCK_METHOD1(Duplicate,
      void(BaseRef));
  MOCK_METHOD0(makeDynamic,
      BaseRef(void));
  MOCK_METHOD2(toStr,
      BaseRef(_List, ));
  MOCK_METHOD0(ObjectClass,
      long(void));
  MOCK_METHOD7(ScanFForVariables,
      void(_AVLList, bool, bool, bool, bool, _AVLListX, long));
  MOCK_METHOD2(ScanFForType,
      void(_SimpleList &, int));
  MOCK_METHOD2(CheckFForDependence,
      bool(, bool));
  MOCK_METHOD1(+,
      operator _Formula(const _Formula));
  MOCK_METHOD1(-,
      operator _Formula(const _Formula));
  MOCK_METHOD1(*,
      operator _Formula(const _Formula));
  MOCK_METHOD1(/,
      operator _Formula(const _Formula));
  MOCK_METHOD1(^,
      operator _Formula(const _Formula));
};
