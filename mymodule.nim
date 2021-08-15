import nimpy

type TestType = ref object of PyNimObjectExperimental
  myField: string

proc setMyField(self: TestType, value: string) {.exportpy.} =
  self.myField = value

proc getMyField(self: TestType): string {.exportpy.} =
  self.myField