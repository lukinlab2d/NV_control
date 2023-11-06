<?xml version='1.0' encoding='UTF-8'?>
<Project Type="Project" LVVersion="13008000">
	<Item Name="My Computer" Type="My Computer">
		<Property Name="server.app.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.control.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="server.tcp.enabled" Type="Bool">false</Property>
		<Property Name="server.tcp.port" Type="Int">0</Property>
		<Property Name="server.tcp.serviceName" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.tcp.serviceName.default" Type="Str">My Computer/VI Server</Property>
		<Property Name="server.vi.callsEnabled" Type="Bool">true</Property>
		<Property Name="server.vi.propertiesEnabled" Type="Bool">true</Property>
		<Property Name="specify.custom.address" Type="Bool">false</Property>
		<Item Name="Device VIs" Type="Folder">
			<Item Name="DeviceClose.vi" Type="VI" URL="../Device VIs/DeviceClose.vi"/>
			<Item Name="DeviceDiscovery.vi" Type="VI" URL="../Device VIs/DeviceDiscovery.vi"/>
			<Item Name="DeviceOpen.vi" Type="VI" URL="../Device VIs/DeviceOpen.vi"/>
			<Item Name="DeviceQuery.vi" Type="VI" URL="../Device VIs/DeviceQuery.vi"/>
			<Item Name="DeviceRead.vi" Type="VI" URL="../Device VIs/DeviceRead.vi"/>
			<Item Name="DeviceWrite.vi" Type="VI" URL="../Device VIs/DeviceWrite.vi"/>
			<Item Name="GetDeviceKeys.vi" Type="VI" URL="../Device VIs/GetDeviceKeys.vi"/>
			<Item Name="OpenFirstValidDevice.vi" Type="VI" URL="../Device VIs/OpenFirstValidDevice.vi"/>
			<Item Name="Shutdown.vi" Type="VI" URL="../Device VIs/Shutdown.vi"/>
		</Item>
		<Item Name="OpenMultipleDevices" Type="Folder">
			<Item Name="OpenMultipleDevices.vi" Type="VI" URL="../OpenMultipleDevices/OpenMultipleDevices.vi"/>
		</Item>
		<Item Name="StepAmplitude" Type="Folder">
			<Item Name="GetStepAmplitudeNegative.vi" Type="VI" URL="../StepAmplitude/GetStepAmplitudeNegative.vi"/>
			<Item Name="GetStepAmplitudePositive.vi" Type="VI" URL="../StepAmplitude/GetStepAmplitudePositive.vi"/>
			<Item Name="SetStepAmplitudeNegative.vi" Type="VI" URL="../StepAmplitude/SetStepAmplitudeNegative.vi"/>
			<Item Name="SetStepAmplitudePositive.vi" Type="VI" URL="../StepAmplitude/SetStepAmplitudePositive.vi"/>
			<Item Name="StepAmplitude.vi" Type="VI" URL="../StepAmplitude/StepAmplitude.vi"/>
		</Item>
		<Item Name="AppendToOutput.vi" Type="VI" URL="../AppendToOutput.vi"/>
		<Item Name="Dependencies" Type="Dependencies">
			<Item Name="CmdLibAgilis.dll" Type="Document" URL="../CmdLibAgilis.dll"/>
			<Item Name="mscorlib" Type="VI" URL="mscorlib">
				<Property Name="NI.PreserveRelativePath" Type="Bool">true</Property>
			</Item>
			<Item Name="VCPIOLib.dll" Type="Document" URL="../VCPIOLib.dll"/>
		</Item>
		<Item Name="Build Specifications" Type="Build"/>
	</Item>
</Project>
